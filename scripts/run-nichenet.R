# run-nichenet.R
# Standard NicheNet — independent per-niche ligand activity analysis
# Per NicheNet v2 protocol paper: single dataset, cell-localization comparison
# Gene set of interest = FindMarkers between the two receiver cell types (run once, split by direction)

run_nichenet <- function(seu_obj) {

    # 1. Load configs
    nn_cfg <- configs$nichenet
    pl_cfg <- pipeline_configs$nichenet
    assay_oi <- configs$paths$DefaultAssay
    lfc_cutoff <- nn_cfg$lfc_cutoff
    expression_pct <- nn_cfg$expression_pct
    top_n_ligands <- nn_cfg$top_n_ligands
    top_n_receptors_per_ligand <- nn_cfg$top_n_receptors_per_ligand
    top_n_targets <- nn_cfg$top_n_targets
    niches <- lapply(nn_cfg$niches, function(n) list(sender = n$sender, receiver = n$receiver))

    # 2. Load NicheNet databases from disk
    org <- nn_cfg$organism
    lr_network <- readRDS(pl_cfg$paths[[org]]$lr_network)
    ligand_target_matrix <- readRDS(pl_cfg$paths[[org]]$ligand_target_matrix)
    lr_network <- lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor)

    # 3. FindMarkers between the two receiver cell types (run once before the per-niche loop)
    #    receiver_ident1 = ident.1 → positive avg_log2FC = upregulated in receiver_ident1
    #    receiver_ident2 = ident.2 → negative avg_log2FC = upregulated in receiver_ident2
    #    logfc.threshold pre-filters to |LFC| >= lfc_cutoff — no low-signal genes make it through
    receiver_ident1 <- niches[[1]]$receiver
    receiver_ident2 <- niches[[2]]$receiver

    message(glue("Computing DE genes: {receiver_ident1} (ident.1) vs {receiver_ident2} (ident.2)"))
    de_between_receivers <- FindMarkers(
        seu_obj,
        ident.1         = receiver_ident1,
        ident.2         = receiver_ident2,
        min.pct         = expression_pct,
        logfc.threshold = lfc_cutoff,
        assay           = assay_oi
    ) %>% filter(p_val_adj < 0.05)

    # Save receiver DE genes: two columns named by cell type, genes fill the rows
    #    Column 1 (receiver_ident1) = genes upregulated in receiver_ident1
    #    Column 2 (receiver_ident2) = genes upregulated in receiver_ident2
    genes_ident1 <- de_between_receivers %>% filter(avg_log2FC > 0) %>% rownames()
    genes_ident2 <- de_between_receivers %>% filter(avg_log2FC < 0) %>% rownames()
    max_len      <- max(length(genes_ident1), length(genes_ident2))
    de_genes_df  <- data.frame(
        a = c(genes_ident1, rep(NA, max_len - length(genes_ident1))),
        b = c(genes_ident2, rep(NA, max_len - length(genes_ident2)))
    )
    colnames(de_genes_df) <- c(receiver_ident1, receiver_ident2)
    write.csv(de_genes_df, glue("{out_base}/niche-net/results/receiver-de-genes.csv"), row.names = FALSE)

    # 4. Per-niche ligand activity analysis
    lapply(names(niches), function(niche_oi) {
        receiver_oi <- niches[[niche_oi]]$receiver
        senders_oi  <- niches[[niche_oi]]$sender

        # 4a. Expressed genes — union across all senders, and for receiver
        expressed_genes_sender <- lapply(senders_oi, function(s) {
            get_expressed_genes(s, seu_obj, pct = expression_pct, assay_oi = assay_oi)
        }) %>% unlist() %>% unique()

        expressed_genes_receiver <- get_expressed_genes(receiver_oi, seu_obj, pct = expression_pct, assay_oi = assay_oi)

        # 4b. Potential ligands (sender-focused): ligands expressed in senders whose receptors are expressed in receiver
        expressed_ligands   <- intersect(lr_network$ligand,   expressed_genes_sender)
        expressed_receptors <- intersect(lr_network$receptor, expressed_genes_receiver)
        potential_ligands   <- lr_network %>%
            filter(ligand %in% expressed_ligands & receptor %in% expressed_receptors) %>%
            pull(ligand) %>% unique()

        # 4c. Gene set of interest: upregulated genes in this niche's receiver vs the other receiver
        #    Since logfc.threshold already enforced |LFC| >= lfc_cutoff, splitting by sign is sufficient
        if (receiver_oi == receiver_ident1) {
            geneset_oi <- de_between_receivers %>% filter(avg_log2FC > 0) %>% rownames()
        } else {
            geneset_oi <- de_between_receivers %>% filter(avg_log2FC < 0) %>% rownames()
        }
        message(glue("{niche_oi} ({receiver_oi}): {length(geneset_oi)} genes in gene set of interest (recommended: 20-1000)"))

        # 4d. Background: receiver expressed genes filtered to ligand_target_matrix gene universe
        #     Per NicheNet v2 vignette: expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
        background_expressed_genes <- expressed_genes_receiver %>%
            .[. %in% rownames(ligand_target_matrix)]

        # 4e. Ligand activity scoring — ranks ligands by how well their known targets predict geneset_oi
        ligand_activities <- predict_ligand_activities(
            geneset                    = geneset_oi,
            background_expressed_genes = background_expressed_genes,
            ligand_target_matrix       = ligand_target_matrix,
            potential_ligands          = potential_ligands
        )
        ligand_activities <- ligand_activities %>%
            arrange(-aupr_corrected) %>%
            mutate(rank = base::rank(-aupr_corrected))

        top_ligands <- ligand_activities %>%
            top_n(top_n_ligands, aupr_corrected) %>%
            pull(test_ligand)

        # 4f. Top target genes per top ligand (genes in geneset_oi predicted downstream of each ligand)
        #     get_weighted_ligand_target_links operates on one ligand at a time — lapply over top_ligands
        active_ligand_target_links_df <- lapply(top_ligands, function(lig) {
            get_weighted_ligand_target_links(
                ligand               = lig,
                geneset              = geneset_oi,
                ligand_target_matrix = ligand_target_matrix,
                n                    = top_n_targets
            )
        }) %>% bind_rows() %>% filter(!is.na(target))

        # 4g. Save per-niche results
        write.csv(ligand_activities,
                  glue("{out_base}/niche-net/results/{niche_oi}-ligand-activity.csv"), row.names = FALSE)
        write.csv(active_ligand_target_links_df,
                  glue("{out_base}/niche-net/results/{niche_oi}-ligand-target.csv"), row.names = FALSE)

        # GO BP on all predicted target genes for this niche (sender-agnostic)
        niche_targets <- active_ligand_target_links_df %>% pull(target) %>% unique()
        run_go_bp(
            genes   = niche_targets,
            label   = niche_oi,
            out_dir = glue("{out_base}/niche-net/go-bp")
        )

        # 4h. Expression tables (ligand, receptor, target) via DotPlot
        #     scale_quantile_adapted scales a numeric vector 0-1 with quantile normalization
        #     mutate_cond caps fraction at expression_pct to prevent outliers from dominating scaling
        features_oi <- union(lr_network$ligand, lr_network$receptor) %>%
            union(active_ligand_target_links_df$target) %>%
            setdiff(NA)

        dotplot <- suppressWarnings(Seurat::DotPlot(
            seu_obj %>% subset(idents = c(senders_oi, receiver_oi) %>% unique()),
            features = features_oi, assay = assay_oi
        ))
        exprs_tbl <- dotplot$data %>% as_tibble() %>%
            dplyr::rename(celltype = id, gene = features.plot, expression = avg.exp,
                   expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
            mutate(fraction = fraction / 100) %>%
            dplyr::select(celltype, gene, expression, expression_scaled, fraction) %>%
            distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))

        exprs_tbl_ligand <- exprs_tbl %>%
            filter(gene %in% lr_network$ligand) %>%
            dplyr::rename(sender = celltype, ligand = gene, ligand_expression = expression,
                   ligand_expression_scaled = expression_scaled, ligand_fraction = fraction) %>%
            mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>%
            mutate(ligand_fraction_adapted = ligand_fraction) %>%
            mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct) %>%
            mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))

        exprs_tbl_receptor <- exprs_tbl %>%
            filter(gene %in% lr_network$receptor) %>%
            dplyr::rename(receiver = celltype, receptor = gene, receptor_expression = expression,
                   receptor_expression_scaled = expression_scaled, receptor_fraction = fraction) %>%
            mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled)) %>%
            mutate(receptor_fraction_adapted = receptor_fraction) %>%
            mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct) %>%
            mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))

        # 4i. Circos plot — top ligands by aupr_corrected, top M receptors per ligand by receptor_fraction
        #     Sender selection is done first in isolation (no receptor columns present) so that
        #     slice_max collapses only the sender dimension. lr_network is then joined against the
        #     resolved sender table, preserving all receptors for the subsequent top_n step.
        best_sender_per_ligand <- exprs_tbl_ligand %>%
            filter(sender %in% senders_oi, ligand %in% top_ligands) %>%
            group_by(ligand) %>%
            slice_max(ligand_fraction, n = 1, with_ties = FALSE) %>%
            ungroup()

        prioritized_tbl_oi <- lr_network %>%
            filter(ligand %in% top_ligands) %>%
            inner_join(best_sender_per_ligand, by = "ligand") %>%
            inner_join(exprs_tbl_receptor %>% filter(receiver == receiver_oi), by = "receptor") %>%
            inner_join(ligand_activities %>% dplyr::select(test_ligand, aupr_corrected),
                       by = c("ligand" = "test_ligand")) %>%
            group_by(ligand) %>%
            top_n(top_n_receptors_per_ligand, receptor_fraction) %>%
            ungroup() %>%
            mutate(
                ligand_receptor      = paste(ligand, receptor, sep = "_"),
                niche                = niche_oi,
                top_niche            = niche_oi,
                prioritization_score = scale_quantile_adapted(aupr_corrected)
            ) %>%
            dplyr::select(niche, sender, receiver, ligand, receptor, ligand_receptor, prioritization_score, top_niche)

        colors_sender <- brewer.pal(
            n    = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% length(),
            name = "Spectral"
        ) %>% magrittr::set_names(prioritized_tbl_oi$sender %>% unique() %>% sort())

        colors_receiver <- c("lavender") %>%
            magrittr::set_names(prioritized_tbl_oi$receiver %>% unique() %>% sort())

        # Open the PDF device BEFORE calling make_circos_lr.
        # circlize draws imperatively to the active device as a side effect — the pdf()
        # call must precede the function or the plot lands in Rplots.pdf instead.
        # separate_legend = TRUE keeps the circos on page 1 clean (base graphics),
        # then we print the legend grob separately on page 2 to avoid cowplot/grid
        # conflicts that arise when replaying a circlize recordedplot in a viewport.
        pdf(glue("{out_base}/niche-net/plots/{niche_oi}-circos.pdf"), width = 12, height = 10)
        circos_output <- make_circos_lr(prioritized_tbl_oi, colors_sender, colors_receiver,
                                        separate_legend = TRUE)
        print(circos_output$p_legend)
        dev.off()
    })
}
