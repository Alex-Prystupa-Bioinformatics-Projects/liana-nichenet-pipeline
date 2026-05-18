# helper-functions.R

# NicheNet helper: conditional mutate — only modifies rows where condition is TRUE
mutate_cond <- function(.data, condition, ...) {
    condition <- eval(substitute(condition), .data, parent.frame())
    .data[condition, ] <- .data[condition, ] %>% mutate(...)
    .data
}

# 0. Create Directories
create.dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = T)
}

create.all.dirs <- function(base_dir) {
   create.dir(base_dir)

   create.dir(glue("{base_dir}/liana"))
   create.dir(glue("{base_dir}/niche-net"))
   create.dir(glue("{base_dir}/liana-nn-combined"))

   create.dir(glue("{base_dir}/liana/plots"))
   create.dir(glue("{base_dir}/niche-net/plots"))
   create.dir(glue("{base_dir}/liana-nn-combined/plots"))

   create.dir(glue("{base_dir}/liana/results"))
   create.dir(glue("{base_dir}/niche-net/results"))
   create.dir(glue("{base_dir}/liana-nn-combined/results"))
   create.dir(glue("{base_dir}/liana-nn-combined/lr-pair-targets"))

   create.dir(glue("{base_dir}/niche-net/go-bp"))
   create.dir(glue("{base_dir}/liana-nn-combined/go-bp"))
}

# 1. Every possible Ligand Sender, Receptor Receiver Combination
get_sender_receiver_combinations <- function(liana_results) {
  cell_types <- unique(c(liana_results$source, liana_results$target))
  cell_type_pairs_df <- expand.grid(cell_types, cell_types) %>% dplyr::rename(Sender = Var1, Receiver = Var2)
  celltype_pairs_list <- lapply(seq_len(nrow(cell_type_pairs_df)), function(i) as.vector(unlist(cell_type_pairs_df[i, ])))

  return(celltype_pairs_list)
}

# 2. LIANA Sender Receiver Filtered Results
get_lr_sr_results <- function(liana_results_df, celltype_pairs_list) {
    
    sr_results_list <- lapply(celltype_pairs_list, function(x) {
        liana_results_df %>% arrange(aggregate_rank) %>%
            filter(aggregate_rank < 0.05) %>%
            filter(source == x[1] & target == x[2]) %>%
    dplyr::rename(ligand=ligand.complex, receptor=receptor.complex)
    })  

    return(sr_results_list)
}

# 3. LIANA Plot Results
# Per Interaction
plot_liana_results_sender_receiver <- function(liana_results, sender, receiver, top_n_interactions=100) {

  # Returns 1 plot
  results <- liana_results %>%
    filter(source == sender & target == receiver) %>%
    dplyr::rename(ligand=ligand.complex, receptor=receptor.complex)

  # filter results to top N interactions
  n <- top_n_interactions
  top_n_results <- results %>%
    arrange(aggregate_rank) %>%
    slice_head(n = n) %>%
    mutate(id = fct_inorder(paste0(ligand, " -> ", receptor)))

  # visualize median rank
  plot <- top_n_results %>%
    ggplot(aes(y = aggregate_rank, x = id)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = configs$combined$aggregate_rank_threshold, color = "red", linetype = "dashed", linewidth = 0.8) +
    xlab("Interaction") + ylab("LIANA's aggregate rank") +
    theme_cowplot() +
    theme(axis.text.x = element_text(size = 8, angle = 60, hjust = 1, vjust = 1))

  
  return(plot)
}

# GO BP enrichment — takes a gene vector, saves a results CSV and dotplot PDF
# label:   used in plot title and output filenames (e.g. "Squamous_niche-T-cells")
# out_dir: directory to write results into (no trailing slash)
run_go_bp <- function(genes, label, out_dir) {
    genes <- unique(genes[!is.na(genes)])
    if (length(genes) < 5) {
        message(glue("Skipping GO BP for {label}: fewer than 5 genes"))
        return(invisible(NULL))
    }

    write.csv(
        data.frame(gene = genes),
        glue("{out_dir}/{label}-go-bp-input-genes.csv"),
        row.names = FALSE
    )

    ego <- enrichGO(
        gene          = genes,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE
    )

    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
        message(glue("No significant GO BP terms for {label}"))
        return(invisible(NULL))
    }

    result_df <- as.data.frame(ego) %>% arrange(p.adjust)

    write.csv(result_df, glue("{out_dir}/{label}-go-bp.csv"), row.names = FALSE)

    p <- result_df %>%
        slice_head(n = 20) %>%
        mutate(Description = fct_reorder(Description, -p.adjust)) %>%
        ggplot(aes(x = Count, y = Description, fill = p.adjust)) +
        geom_bar(stat = "identity") +
        scale_fill_gradient(low = "#d73027", high = "#f7f7f7", name = "Adj. p-value") +
        xlab("Gene count") +
        ylab(NULL) +
        ggtitle(glue("GO BP: {label}")) +
        theme_cowplot() +
        theme(axis.text.y = element_text(size = 8))

    pdf(glue("{out_dir}/{label}-go-bp.pdf"), width = 10, height = 8)
    print(p)
    dev.off()
}

# Plot All Results
plot_liana_all <- function(celltype_pairs_list, liana_results) {

  # Get Plots
  plot_list <- lapply(celltype_pairs_list, function(x) {
    plot_liana_results_sender_receiver(liana_results, sender = x[1], receiver = x[2])
  })

  plot_names <- lapply(celltype_pairs_list, function(x) {
    paste0(x[1], "_", x[2]) %>% unlist()
  })

  names(plot_list) <- plot_names
  
  return(plot_list)
}
