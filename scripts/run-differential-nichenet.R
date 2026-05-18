# run-differential-nichenet.R

run_differential_nichenet <- function(seu_obj) {

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

    # Prioritization weights — see README.md for explanation of each weight
    prioritizing_weights <- c(
        scaled_ligand_score                        = 5,    # niche-specific DE of ligand (most important)
        scaled_ligand_expression_scaled            = 1,    # scaled ligand expression across cell types
        ligand_fraction                            = 1,    # fraction of sender cells expressing the ligand
        scaled_ligand_score_spatial                = 0,    # spatial DE of ligand (no spatial data — zeroed)
        scaled_receptor_score                      = 0.5,  # niche-specific DE of receptor
        scaled_receptor_expression_scaled          = 0.5,  # scaled receptor expression across cell types
        receptor_fraction                          = 1,    # fraction of receiver cells expressing the receptor
        ligand_scaled_receptor_expression_fraction = 1,    # receptor expression strength within ligand family
        scaled_receptor_score_spatial              = 0,    # spatial DE of receptor (no spatial data — zeroed)
        scaled_activity                            = 0,    # absolute ligand activity (not recommended raw)
        scaled_activity_normalized                 = 1     # normalized ligand activity
    )

    # 2. Load NicheNet databases from disk
    org <- nn_cfg$organism
    lr_network <- readRDS(pl_cfg$paths[[org]]$lr_network)
    ligand_target_matrix <- readRDS(pl_cfg$paths[[org]]$ligand_target_matrix)
    lr_network <- lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor)

    # 3. Differential expression — senders and receivers across niches
    DE_sender <- calculate_niche_de(seurat_obj = seu_obj, niches = niches, type = "sender", assay_oi = assay_oi)
    DE_sender_processed <- process_niche_de(DE_table = DE_sender, niches = niches, expression_pct = expression_pct, type = "sender")

    DE_receiver <- calculate_niche_de(seurat_obj = seu_obj, niches = niches, type = "receiver", assay_oi = assay_oi)
    DE_receiver_processed <- process_niche_de(DE_table = DE_receiver, niches = niches, expression_pct = expression_pct, type = "receiver")

    DE_sender_receiver <- combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = "min_lfc")

    # 4. Non-spatial DE — zeroed placeholder tables (no spatial data; spatial weights = 0)
    spatial_info <- tibble(celltype_region_oi = NA, celltype_other_region = NA, niche = NA, celltype_type = NA)
    sender_spatial_DE_processed <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
    sender_spatial_DE_processed <- sender_spatial_DE_processed %>%
        mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))

    receiver_spatial_DE_processed <- get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
    receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>%
        mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))

    # 5. Target DE per niche
    DE_receiver_targets <- calculate_niche_de_targets(seurat_obj = seu_obj, niches = niches,
                                                       expression_pct = expression_pct, lfc_cutoff = lfc_cutoff,
                                                       assay_oi = assay_oi)
    DE_receiver_processed_targets <- process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, niches = niches,
                                                                 expression_pct = expression_pct,
                                                                 specificity_score = "min_lfc")

    # 6. Build geneset of interest per niche
    #    Geneset = receiver genes upregulated above lfc_cutoff in this niche vs all others
    #    Vignette recommends geneset size 20-1000 genes; adjust lfc_cutoff in config.yml if outside range
    niche_geneset_list <- lapply(names(niches), function(niche_oi) {
        receiver_oi <- niches[[niche_oi]]$receiver

        geneset <- DE_receiver_processed_targets %>%
            filter(receiver == receiver_oi & niche == niche_oi & target_score > lfc_cutoff) %>%
            pull(target) %>% unique()

        background <- rownames(seu_obj)[Matrix::rowSums(GetAssayData(seu_obj, assay = assay_oi)) > 0]

        list(receiver = receiver_oi, geneset = geneset, background = background)
    })
    names(niche_geneset_list) <- names(niches)

    # 7. Ligand activity calculation
    ligand_activities_targets <- get_ligand_activities_targets(
        niche_geneset_list = niche_geneset_list,
        ligand_target_matrix = ligand_target_matrix,
        top_n_target = top_n_targets
    )

    # 8. Expression tables (ligand, receptor, target) via DotPlot
    #    scale_quantile_adapted scales a numeric vector 0-1 with quantile normalization
    #    mutate_cond caps fraction at expression_pct to prevent outliers from dominating scaling
    features_oi <- union(lr_network$ligand, lr_network$receptor) %>%
        union(ligand_activities_targets$target) %>%
        setdiff(NA)

    dotplot <- suppressWarnings(Seurat::DotPlot(
        seu_obj %>% subset(idents = niches %>% unlist() %>% unique()),
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

    exprs_tbl_target <- exprs_tbl %>%
        filter(gene %in% ligand_activities_targets$target) %>%
        dplyr::rename(receiver = celltype, target = gene, target_expression = expression,
               target_expression_scaled = expression_scaled, target_fraction = fraction)

    # Build ligand_scaled_receptor_expression_fraction_df
    #    Ranks receptors by expression and fraction within each ligand family, then combines into a 0-1 score
    exprs_sender_receiver <- lr_network %>%
        inner_join(exprs_tbl_ligand, by = c("ligand")) %>%
        inner_join(exprs_tbl_receptor, by = c("receptor")) %>%
        inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))

    ligand_scaled_receptor_expression_fraction_df <- exprs_sender_receiver %>%
        group_by(ligand, receiver) %>%
        mutate(rank_receptor_expression = dense_rank(receptor_expression),
               rank_receptor_fraction   = dense_rank(receptor_fraction)) %>%
        mutate(ligand_scaled_receptor_expression_fraction = 0.5 * (
            (rank_receptor_fraction  / max(rank_receptor_fraction)) +
            (rank_receptor_expression / max(rank_receptor_expression))
        )) %>%
        distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction) %>%
        distinct() %>% ungroup()

    # 9. Prioritization — combines all scores into final ranked LR and LT tables
    output <- list(
        DE_sender_receiver                          = DE_sender_receiver,
        ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df,
        sender_spatial_DE_processed                 = sender_spatial_DE_processed,
        receiver_spatial_DE_processed               = receiver_spatial_DE_processed,
        ligand_activities_targets                   = ligand_activities_targets,
        DE_receiver_processed_targets               = DE_receiver_processed_targets,
        exprs_tbl_ligand                            = exprs_tbl_ligand,
        exprs_tbl_receptor                          = exprs_tbl_receptor,
        exprs_tbl_target                            = exprs_tbl_target
    )

    prioritization_tables <- get_prioritization_tables(output, prioritizing_weights)

    # Save Results
    write.csv(prioritization_tables$prioritization_tbl_ligand_receptor,
              glue("{out_base}/niche-net/results/nichenet-prioritization-lr.csv"))
    write.csv(prioritization_tables$prioritization_tbl_ligand_target,
              glue("{out_base}/niche-net/results/nichenet-prioritization-lt.csv"))

    # 10. Circos plots — one per niche, colors generated dynamically via brewer.pal
    top_ligand_receptor_niche_df <- prioritization_tables$prioritization_tbl_ligand_receptor %>%
        dplyr::select(niche, sender, receiver, ligand, receptor, prioritization_score) %>%
        group_by(ligand, receptor) %>% top_n(1, prioritization_score) %>%
        ungroup() %>% dplyr::select(ligand, receptor, niche) %>%
        dplyr::rename(top_niche = niche)

    for (niche_oi in names(niches)) {
        receiver_oi <- niches[[niche_oi]]$receiver

        # Top N ligands for this receiver
        filtered_ligands <- prioritization_tables$prioritization_tbl_ligand_receptor %>%
            filter(receiver == receiver_oi) %>%
            top_n(top_n_ligands, prioritization_score) %>%
            pull(ligand) %>% unique()

        # Top M receptors per ligand, joined with top niche info
        prioritized_tbl_oi <- prioritization_tables$prioritization_tbl_ligand_receptor %>%
            filter(ligand %in% filtered_ligands) %>%
            dplyr::select(niche, sender, receiver, ligand, receptor, ligand_receptor, prioritization_score) %>%
            distinct() %>%
            inner_join(top_ligand_receptor_niche_df) %>%
            group_by(ligand) %>%
            filter(receiver == receiver_oi) %>%
            top_n(top_n_receptors_per_ligand, prioritization_score) %>%
            ungroup()

        colors_sender <- brewer.pal(
            n = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% length(),
            name = "Spectral"
        ) %>% magrittr::set_names(prioritized_tbl_oi$sender %>% unique() %>% sort())

        colors_receiver <- c("lavender") %>% magrittr::set_names(prioritized_tbl_oi$receiver %>% unique() %>% sort())

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
    }
}
