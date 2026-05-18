# run-liana.R

run_liana <- function(seu_obj) {

    # Full Results
    liana_results <- liana_wrap(seu_obj) %>% liana_aggregate()

    # Sender -> Receiver Filtered Specific Results
    sender_receiver_combinations <- get_sender_receiver_combinations(liana_results)
    lr_sr_results <- get_lr_sr_results(liana_results, sender_receiver_combinations)

    # Bar Plots
    liana_barplots <- plot_liana_all(sender_receiver_combinations, liana_results)

    # Save
    # 1. Save Results
    write.csv(liana_results, glue("{out_base}/liana/results/liana-lr-full-results.csv"))

    for (sr in sender_receiver_combinations) {
        results <- liana_results %>%
            filter(source == sr[1] & target == sr[2]) %>%
            dplyr::rename(ligand=ligand.complex, receptor=receptor.complex)

        plot_name <- paste0(sr[1], "_", sr[2])

        write.csv(results, glue("{out_base}/liana/results/liana-lr-{plot_name}-results.csv"))
    }


    # 2. Save Plots
    for (sr in names(liana_barplots)) {
        plot <- liana_barplots[[sr]]
        pdf(glue("{out_base}/liana/plots/{sr}.pdf"), height = 8, width = 12)
        print(plot)
        dev.off()
    }
}
