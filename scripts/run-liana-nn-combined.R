# run-liana-nn-combined.R
# Implements the LIANA + NicheNet combined approach from the liana vignette:
# https://saezlab.github.io/liana/articles/liana_nichenet.html
#
# Per niche, per sender:
#   1. Filter LIANA to this sender→receiver pair only
#   2. Extract unique LIANA ligands for this sender, filter to those in ligand_target_matrix
#   3. Re-run predict_ligand_activities using only those sender-specific LIANA ligands
#   4. Take top N by pearson, inner-join back to LIANA
#   5. Barplot (NicheNet activity) + heatmap (LIANA LR pairs) — one PDF per sender→receiver
#   6. GO BP on target genes for this sender
#
# After per-sender loop, per niche:
#   7. Circos using LIANA aggregate_rank as prioritization_score — per-(sender,ligand,receptor)
#      scoring gives true sender differentiation in the circos

run_liana_nn_combined <- function(seu_obj) {

    # 1. Load configs
    combined_cfg   <- configs$combined
    niches         <- configs$nichenet$niches
    top_n_liana    <- combined_cfg$top_n_liana_interactions
    rank_thresh    <- combined_cfg$aggregate_rank_threshold
    top_n_nn       <- combined_cfg$top_n_nn_ligands
    top_n_ligands  <- configs$nichenet$top_n_ligands
    top_n_targets  <- configs$nichenet$top_n_targets
    top_n_receptors_per_ligand <- configs$nichenet$top_n_receptors_per_ligand
    expression_pct <- configs$nichenet$expression_pct
    assay_oi       <- configs$paths$DefaultAssay

    # 2. Load LIANA full results once (outside loop)
    liana_full <- read.csv(glue("{out_base}/liana/results/liana-lr-full-results.csv")) %>%
        dplyr::rename(ligand = ligand.complex, receptor = receptor.complex)

    # 3. Load NicheNet ligand_target_matrix
    org <- configs$nichenet$organism
    ligand_target_matrix <- readRDS(pipeline_configs$nichenet$paths[[org]]$ligand_target_matrix)

    # 4. Load pre-computed receiver DE genes (gene sets of interest, saved by run-nichenet.R)
    de_genes <- read.csv(glue("{out_base}/niche-net/results/receiver-de-genes.csv"))

    # 5. Per-niche loop
    lapply(names(niches), function(niche_oi) {

        receiver_oi <- niches[[niche_oi]]$receiver
        senders_oi  <- niches[[niche_oi]]$sender

        message(glue("Combined analysis: {niche_oi} (receiver: {receiver_oi})"))

        # 5a. Gene set of interest for this receiver (from saved CSV)
        geneset_oi <- de_genes[[receiver_oi]] %>% na.omit() %>% as.character()
        message(glue("{niche_oi}: {length(geneset_oi)} genes in gene set of interest"))

        # 5b. Background genes: receiver expressed genes filtered to ligand_target_matrix gene universe
        #     Per NicheNet v2 vignette: expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
        background_genes <- get_expressed_genes(receiver_oi, seu_obj,
                                                pct = expression_pct, assay_oi = assay_oi) %>%
            .[. %in% rownames(ligand_target_matrix)]

        # 5b-diag. Diagnostic: gene set vs background, with LTM membership
        effective_geneset <- intersect(geneset_oi, rownames(ligand_target_matrix))
        message(glue("{niche_oi}: {length(geneset_oi)} DE genes -> {length(effective_geneset)} effective in LTM ({round(100*length(effective_geneset)/max(length(geneset_oi),1),1)}%)"))
        message(glue("{niche_oi}: background = {length(background_genes)} genes"))

        max_diag <- max(length(geneset_oi), length(background_genes))
        diag_df <- data.frame(
            de_gene         = c(geneset_oi,
                                rep(NA, max_diag - length(geneset_oi))),
            de_in_ltm       = c(geneset_oi %in% rownames(ligand_target_matrix),
                                rep(NA, max_diag - length(geneset_oi))),
            effective_gene  = c(effective_geneset,
                                rep(NA, max_diag - length(effective_geneset))),
            background_gene = c(background_genes,
                                rep(NA, max_diag - length(background_genes)))
        )
        write.csv(diag_df,
                  glue("{out_base}/liana-nn-combined/results/{niche_oi}-geneset-background-diagnostic.csv"),
                  row.names = FALSE)

        # 5c. Per-sender loop — independent barplot + heatmap per sender→receiver pair
        # Collect nichenet_activities per sender so 5d can use aupr_corrected in the circos
        sender_activities <- list()

        lapply(senders_oi, function(sender_oi) {

            message(glue("  Sender: {sender_oi} -> {receiver_oi}"))

            # 5c-i. Filter LIANA to this sender→receiver pair only
            liana_sender <- liana_full %>%
                filter(source == sender_oi & target == receiver_oi) %>%
                filter(aggregate_rank < rank_thresh) %>%
                arrange(aggregate_rank) %>%
                slice_head(n = top_n_liana)

            # 5c-ii. Extract ligands for this sender, filter to those in ligand_target_matrix
            ligands <- unique(liana_sender$ligand)
            ligands <- ligands[ligands %in% colnames(ligand_target_matrix)]

            if (length(ligands) == 0) {
                message(glue("  WARNING: No LIANA ligands in ligand_target_matrix for {sender_oi} -> {receiver_oi} — skipping."))
                return(NULL)
            }

            message(glue("  {sender_oi}: {length(ligands)} LIANA ligands in ligand_target_matrix"))

            # 5c-iii. NicheNet activity — scored on this sender's ligands only
            nichenet_activities <- predict_ligand_activities(
                geneset                    = geneset_oi,
                background_expressed_genes = background_genes,
                ligand_target_matrix       = ligand_target_matrix,
                potential_ligands          = ligands
            ) %>% arrange(desc(pearson))

            # Save for use in 5d circos
            sender_activities[[sender_oi]] <<- nichenet_activities

            # 5c-iv. Top N ligands by pearson, join back to LIANA
            top_nn_ligands <- nichenet_activities %>%
                slice_head(n = top_n_nn) %>%
                pull(test_ligand)

            liana_nn_joined <- liana_sender %>%
                filter(ligand %in% top_nn_ligands) %>%
                inner_join(
                    nichenet_activities %>% dplyr::select(test_ligand, auroc, aupr_corrected, pearson),
                    by = c("ligand" = "test_ligand")
                ) %>%
                arrange(pearson) %>%
                mutate(ligand = fct_inorder(ligand))

            ligand_order <- levels(liana_nn_joined$ligand)

            if (nrow(liana_nn_joined) == 0) {
                message(glue("  WARNING: Empty join for {sender_oi} -> {receiver_oi} — skipping."))
                return(NULL)
            }

            # 5c-v. Compute target genes for the top ligands
            top_ligands_chr <- unique(liana_nn_joined$ligand) %>% as.character()
            nn_targets <- lapply(top_ligands_chr, function(lig) {
                get_weighted_ligand_target_links(
                    ligand               = lig,
                    geneset              = geneset_oi,
                    ligand_target_matrix = ligand_target_matrix,
                    n                    = top_n_targets
                )
            }) %>% bind_rows() %>%
                filter(!is.na(target)) %>%
                dplyr::rename(target_gene = target, target_weight = weight)

            liana_nn_joined <- liana_nn_joined %>%
                left_join(nn_targets, by = "ligand", relationship = "many-to-many")

            # 5c-vi. Save per-sender CSV
            # Sanitize sender name for use in file paths (replace special chars with _)
            sender_safe <- gsub("[^A-Za-z0-9_]", "_", sender_oi)
            write.csv(liana_nn_joined,
                      glue("{out_base}/liana-nn-combined/results/{niche_oi}-{sender_safe}-liana-nn-combined.csv"),
                      row.names = FALSE)

            # 5c-vi-b. Save per-LR-pair target gene CSV (wide format)
            # Columns = LR pair (ligand_receptor), rows = predicted target genes
            lr_pair_wide <- liana_nn_joined %>%
                filter(!is.na(target_gene)) %>%
                mutate(lr_pair = paste(ligand, receptor, sep = "_")) %>%
                dplyr::select(lr_pair, target_gene, target_weight) %>%
                distinct() %>%
                arrange(lr_pair, desc(target_weight)) %>%
                group_by(lr_pair) %>%
                mutate(row = row_number()) %>%
                ungroup() %>%
                dplyr::select(-target_weight) %>%
                tidyr::pivot_wider(names_from = lr_pair, values_from = target_gene) %>%
                dplyr::select(-row)

            write.csv(lr_pair_wide,
                      glue("{out_base}/liana-nn-combined/lr-pair-targets/{niche_oi}-{sender_safe}-lr-pair-targets.csv"),
                      row.names = FALSE)

            # 5c-vii. Plot A: NicheNet barplot
            p_bar <- nichenet_activities %>%
                filter(test_ligand %in% top_nn_ligands) %>%
                mutate(ligand = factor(test_ligand, levels = ligand_order)) %>%
                ggplot(aes(x = pearson, y = ligand)) +
                geom_bar(stat = "identity", fill = "steelblue") +
                xlab("NicheNet ligand activity\n(Pearson's score)") +
                ylab("Ligand") +
                theme_cowplot() +
                ggtitle(glue("{niche_oi} | {sender_oi}\nNicheNet Activity"))

            # 5c-viii. Plot B: LIANA heatmap — single sender, no faceting
            p_heat <- liana_nn_joined %>%
                mutate(
                    ligand   = factor(ligand,   levels = ligand_order),
                    receptor = factor(receptor, levels = sort(unique(receptor)))
                ) %>%
                ggplot(aes(x = receptor, y = ligand, fill = aggregate_rank)) +
                geom_tile(color = "white", linewidth = 0.3) +
                scale_fill_gradient(
                    low  = "#d73027",
                    high = "#f7f7f7",
                    name = "Aggregate rank\n(lower = more significant)"
                ) +
                theme_cowplot() +
                theme(
                    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7),
                    axis.text.y = element_text(size = 7)
                ) +
                xlab("Receptor") +
                ylab("Ligand") +
                ggtitle(glue("{niche_oi} | {sender_oi}\nSignificant LIANA interactions"))

            # 5c-ix. Combine and save per-sender PDF
            p_combined <- plot_grid(
                p_bar, p_heat,
                nrow       = 1,
                rel_widths = c(1, 2),
                align      = "h",
                axis       = "tb"
            )

            pdf(glue("{out_base}/liana-nn-combined/plots/{niche_oi}-{sender_safe}-liana-nn-combined.pdf"),
                width  = 14,
                height = max(8, n_distinct(liana_nn_joined$ligand) * 0.4))
            print(p_combined)
            dev.off()

            # 5c-x. GO BP on target genes for this sender
            # Only use target genes from ligands with positive pearson — negative-pearson
            # ligands have no predictive power for receiver expression and should not
            # contribute to enrichment analysis.
            positive_pearson_ligands <- nichenet_activities %>%
                filter(pearson > 0, test_ligand %in% top_ligands_chr) %>%
                pull(test_ligand)

            sender_targets <- nn_targets %>%
                filter(ligand %in% positive_pearson_ligands) %>%
                pull(target_gene) %>% unique()
            run_go_bp(
                genes   = sender_targets,
                label   = glue("{niche_oi}-{sender_safe}"),
                out_dir = glue("{out_base}/liana-nn-combined/go-bp")
            )
        })

        # 5d. Niche-level circos — all senders combined, aupr_corrected as prioritization score
        #     aupr_corrected is joined per sender from sender_activities (collected in 5c)
        liana_niche_all <- lapply(senders_oi, function(sender_oi) {
            liana_full %>%
                filter(source == sender_oi & target == receiver_oi) %>%
                filter(aggregate_rank < rank_thresh) %>%
                arrange(aggregate_rank) %>%
                slice_head(n = top_n_liana) %>%
                left_join(
                    sender_activities[[sender_oi]] %>%
                        dplyr::select(test_ligand, aupr_corrected),
                    by = c("ligand" = "test_ligand")
                )
        }) %>% bind_rows()

        if (nrow(liana_niche_all) == 0) {
            message(glue("WARNING: No LIANA interactions for {niche_oi} circos — skipping."))
            return(NULL)
        }

        # Cap to top_n_ligands by best aggregate_rank — LIANA determines which interactions to show
        top_ligands_circos <- liana_niche_all %>%
            group_by(ligand) %>%
            summarise(best_rank = min(aggregate_rank), .groups = "drop") %>%
            slice_min(best_rank, n = top_n_ligands, with_ties = FALSE) %>%
            pull(ligand)
        liana_niche_all <- liana_niche_all %>%
            filter(ligand %in% top_ligands_circos) %>%
            filter(!is.na(aupr_corrected))

        # Shift aupr_corrected so the most-negative value becomes 0.05.
        # This makes arc sizes proportional to NicheNet activity rather than
        # inflating arcs for negative-aupr ligands with no visible chord.
        # No shift applied when all values are already positive.
        min_aupr <- min(liana_niche_all$aupr_corrected, na.rm = TRUE)
        shift    <- if (min_aupr < 0) abs(min_aupr) + 0.05 else 0
        liana_niche_all <- liana_niche_all %>%
            mutate(aupr_corrected_orig = aupr_corrected,
                   aupr_corrected      = aupr_corrected + shift)

        prioritized_tbl_circos <- liana_niche_all %>%
            dplyr::select(source, target, ligand, receptor, aggregate_rank, aupr_corrected, aupr_corrected_orig) %>%
            dplyr::rename(sender = source, receiver = target) %>%
            mutate(ligand_receptor = paste(ligand, receptor, sep = "_"),
                   niche           = niche_oi,
                   top_niche       = niche_oi) %>%
            group_by(ligand) %>%
            top_n(top_n_receptors_per_ligand, -aggregate_rank) %>%    # receptor selection: same as before
            ungroup() %>%
            mutate(prioritization_score = aupr_corrected)             # line width: shifted aupr_corrected

        # Chords where aupr_corrected was originally negative are fully transparent
        # (no visible arrow), while originally-positive chords are fully opaque.
        transparency_vec <- ifelse(prioritized_tbl_circos$aupr_corrected_orig < 0, 1, 0)

        prioritized_tbl_circos <- prioritized_tbl_circos %>%
            dplyr::select(niche, sender, receiver, ligand, receptor,
                          ligand_receptor, prioritization_score, top_niche)

        if (nrow(prioritized_tbl_circos) > 0) {
            colors_sender_circos <- brewer.pal(
                n    = prioritized_tbl_circos$sender %>% unique() %>% sort() %>% length(),
                name = "Spectral"
            ) %>% magrittr::set_names(prioritized_tbl_circos$sender %>% unique() %>% sort())

            colors_receiver_circos <- c("lavender") %>%
                magrittr::set_names(prioritized_tbl_circos$receiver %>% unique() %>% sort())

            pdf(glue("{out_base}/liana-nn-combined/plots/{niche_oi}-liana-nn-combined-circos.pdf"),
                width = 12, height = 10)
            circos_output <- make_circos_lr(prioritized_tbl_circos,
                                            colors_sender_circos,
                                            colors_receiver_circos,
                                            transparency = transparency_vec,
                                            separate_legend = TRUE)
            print(circos_output$p_legend)
            dev.off()
        }
    })
}
