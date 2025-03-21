#model = rma.uni or mv model
#mod is moderator (can by 1 if intercept only)
#es is the name of effect size - for the plot
#group is the grouping variable for the plot
#terms in case you want to specify your own labels and not default output. 
#pb if pub bias = components you'd liek to hold constant,
#it will automatically correct the estimates to when pub bias components are zero
orchaRd_table <- function(model, mod, es, group, terms = NULL, pb = NULL) {
  
  require(orchaRd)
  require(tidyverse)
  require(broom)
  require(patchwork)
  require(gt)
  
  model_summary <- orchaRd::mod_results(model, mod = mod, group = group)
  
  # orchard plot
  if(!is.null(pb)){
  correct <-setNames(as.list(rep(0, length(pb))), pb)
  model_summary <- orchaRd::mod_results(model, mod = mod, group = group,
                                        at = correct)
  }
  
  orchard <- orchaRd::orchard_plot(model_summary, xlab = es)
  
  print(model_summary)
  
  
  # table for summary
  main_table <- broom::tidy(model, conf.int = TRUE)
  
  if(is.null(terms)){
    terms = main_table$term
  }
  
  main_table <- main_table %>%
    dplyr::mutate(
      Predictor = terms,
      Estimate = sprintf("%.3f", estimate),
      SE = sprintf("%.3f", std.error),
      CI.Low = sprintf("%.3f", conf.low),
      CI.High = sprintf("%.3f", conf.high),
      Stat = sprintf("%.3f", statistic),
      p = sprintf("%.3f", p.value)
    ) %>%
    dplyr::mutate(p = case_when(
      p < 0.001 ~ "<0.001",
      p >= 0.001 ~ p
    )) %>%
    select(
      Predictor,
      Estimate,
      SE,
      CI.Low,
      CI.High,
      Stat,
      p
    )
  
  layout <- "A#
             BB
             BB"
  
  combined_plot <-
    wrap_table(
      main_table |>
        gt() |>
        tab_style(
          style = list(cell_text(weight = "bold")),
          locations = list(cells_body(rows = p < 0.05))
        ) |>
        tab_style(
          style = list(cell_text(style = "italic")),
          locations = list(cells_column_labels())
        ) |>
        cols_label(
          Predictor = "Predictor",
          Estimate = "Estimate",
          SE = "SE",
          CI.Low = "CI Low",
          CI.High = "CI High",
          Stat = "Statistic",
          p = "p"
        ) |>
        tab_style(
          style = cell_borders(
            sides = "bottom",
            color = "black",
            weight = px(6) # Thickness of 3 pixels
          ),
          locations = cells_column_labels()
        ) |>
        cols_align("left") |>
        opt_table_font(font = "Arial") |>
        tab_options(
          table.width = px(750),
          table.border.bottom.width = px(0),
          source_notes.border.bottom.width = px(0)
        ) |>
        tab_source_note(
          source_note = paste(
            "Test for Residual Heterogeneity: QE(df = ",
            round(model$QEdf, 3), ") = ", round(model$QE, 3),
            ", p = ", ifelse(model$QEp < 0.001, "<0.001", round(model$QEp, 3)), 
            "\nTest for Moderators: F(df1 = ",
            model$QMdf[1], ", df2 = ", model$QMdf[2],
            ") = ", round(model$QM, 3),
            ", p = ", ifelse(model$QMp < 0.001, "<0.001", round(model$QMp, 3)),
            sep = ""
          )
        )
    ) / orchard +
    plot_layout(heights = c(2, 6), design = layout)
  
  return(combined_plot)
}
