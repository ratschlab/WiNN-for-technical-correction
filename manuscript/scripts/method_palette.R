# One colour per method, shared by every main and supplementary figure.
#
# Each method keeps the same colour in every panel, so a method can be recognised
# without reading its label. All colours are saturated enough to read on white,
# including as text labels. WiNN is red; its two supplementary configurations use
# a dark red and a magenta. Panels that show only WiNN (gate policies, alignment
# axes, workflow stages) draw everything in the WiNN colour and separate the
# conditions by shape or line type, so colour always means "method".
#
# Source from a figure script with:
#   source(file.path(dirname(script_path), "method_palette.R"))

palette_method_ids <- c(
  "raw", "combat", "qc_rlsc", "qc_rfsc", "tiger", "serrf",
  "winn_auto_qc", "winn_auto_batch_qc", "winn_fixed_no_qc"
)

palette_method_labels <- c(
  raw = "Raw", combat = "ComBat", qc_rlsc = "QC-RLSC", qc_rfsc = "QC-RFSC",
  tiger = "TIGER", serrf = "SERRF",
  winn_auto_qc = "WiNN (QC-tuned)",
  winn_auto_batch_qc = "WiNN (inferred segments)",
  winn_fixed_no_qc = "WiNN"
)

palette_method_colours <- c(
  "Raw" = "#6E6E6E",
  "ComBat" = "#0072B2",
  "QC-RLSC" = "#009473",
  "QC-RFSC" = "#6A3D9A",
  "TIGER" = "#D27500",
  "SERRF" = "#7F4F24",
  "WiNN (QC-tuned)" = "#8B0000",
  "WiNN (inferred segments)" = "#E7298A",
  "WiNN" = "#D7191C"
)

# Colours keyed by internal method identifier.
palette_colours_by_id <- setNames(
  unname(palette_method_colours[palette_method_labels[palette_method_ids]]),
  palette_method_ids
)

winn_colour <- palette_method_colours[["WiNN"]]

# Outline and error-bar ink shared by all panels.
palette_ink <- "#172033"

# Seed- or split-level points are drawn in the method colour at this opacity;
# summary points (means, point estimates) are drawn fully opaque.
palette_point_alpha <- 0.55

# Shapes for conditions within WiNN-only panels.
gate_shapes <- c(
  "Selective / selective" = 16, "Forced drift" = 17,
  "Forced batch" = 15, "Forced drift + batch" = 18
)
alignment_shapes <- c("Batch" = 16, "Run order" = 17, "Both" = 15)
alignment_linetypes <- c("Batch" = "solid", "Run order" = "dashed", "Both" = "dotdash")
