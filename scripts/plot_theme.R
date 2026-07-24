# Shared plotting conventions for benchmark figures.

winn_method_order <- function() {
  c(
    "Raw", "ComBat", "QC-RLSC", "QC-RFSC", "TIGER", "SERRF",
    "WINN auto (QC)", "WINN auto-batch (QC)", "WINN default (no QC)"
  )
}

winn_method_palette <- function() {
  c(
    "Raw" = "#5F6368",
    "ComBat" = "#0072B2",
    "QC-RLSC" = "#009E73",
    "QC-RFSC" = "#56B4E9",
    "TIGER" = "#E69F00",
    "SERRF" = "#D55E00",
    "WINN auto (QC)" = "#CC79A7",
    "WINN auto-batch (QC)" = "#6F4E7C",
    "WINN default (no QC)" = "#B79F00"
  )
}

winn_dataset_palette <- function() {
  c(
    "Simulation" = "#4477AA",
    "MTBLS79" = "#228833",
    "BatchCorr Set 1" = "#CC6677",
    "Sacurine" = "#AA3377",
    "WaveICA" = "#66CCEE"
  )
}

winn_ablation_gate_palette <- function() {
  c(
    "Selective/Selective" = "#4477AA",
    "All/Selective" = "#EE6677",
    "Selective/All" = "#228833",
    "All/All" = "#AA3377"
  )
}

winn_selection_palette <- function() {
  c(
    "Drift-corrected feature-batch profiles" = "#4477AA",
    "Features corrected for drift" = "#66CCEE",
    "Features corrected for residual batch" = "#228833",
    "Samples altered by PQN shrinkage" = "#CCBB44"
  )
}

winn_theme_publication <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size, base_family = "Helvetica") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(
        colour = "#E2E5E9", linewidth = 0.35, linetype = "dotted"
      ),
      axis.text = ggplot2::element_text(colour = "#343A40"),
      axis.title = ggplot2::element_text(colour = "#202428", face = "bold"),
      strip.background = ggplot2::element_rect(fill = "#EEF1F4", colour = NA),
      strip.text = ggplot2::element_text(colour = "#202428", face = "bold"),
      plot.title = ggplot2::element_text(colour = "#202428", face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "#59636E"),
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold")
    )
}
