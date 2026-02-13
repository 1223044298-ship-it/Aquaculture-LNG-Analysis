library(tidyverse)
library(cowplot)

# =========================
# Output directory (relative path)
# =========================
output_dir <- file.path("results", "fig2_radial")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# Palette & Order
# =========================
aa_colors <- c(
  "Ile" = "#41b6c4", "Leu" = "#225ea8", "Val" = "#a1dab4",
  "AAA" = "#fecc5c", "Trp" = "#fd8d3c",
  "His" = "#9e9ac8", "Lys" = "#54278f",
  "Thr" = "#f03b20", "SAA" = "#bdbdbd"
)

aa_order <- c("His","Ile","Leu","Lys","SAA","AAA","Thr","Trp","Val")

# =========================
# Data
# =========================
df <- tribble(
  ~AA,  ~Open_Flow, ~Dense_Array, ~Semi_Enclosed,
  "His", 336.18,  840.22, 2676.28,
  "Ile", 561.55,  1403.50, 4470.43,
  "Leu", 974.80,  2436.36, 7760.31,
  "Lys", 1106.99, 2766.73, 8812.63,
  "SAA", 491.05,  1227.29, 3909.18,
  "AAA", 1015.52, 2538.12, 8084.45,
  "Thr", 417.15,  1042.60, 3320.91,
  "Trp", 131.35,  328.29, 1045.68,
  "Val", 622.85,  1556.71, 4958.46
)

rec_df <- tribble(
  ~AA, ~rec,
  "His", 535.18, "Ile", 1060.58, "Leu", 2070.64,
  "Lys", 1601.10, "SAA", 799.73, "AAA", 1340.24,
  "Thr", 804.21, "Trp", 214.48, "Val", 1378.94
) %>%
  mutate(AA = factor(AA, levels = aa_order)) %>%
  arrange(AA) %>%
  mutate(id = row_number())

df_long <- df %>%
  pivot_longer(cols = c(Open_Flow, Dense_Array, Semi_Enclosed),
               names_to = "scenario", values_to = "value") %>%
  mutate(AA = factor(AA, levels = aa_order)) %>%
  arrange(scenario, AA) %>%
  group_by(scenario) %>%
  mutate(id = row_number()) %>%
  ungroup()

# =========================
# Plot settings
# =========================
ymax <- max(c(df_long$value, rec_df$rec), na.rm = TRUE)
ymin_inner <- -0.5 * ymax
ring_y <- -0.02 * ymax
breaks_y <- seq(2000, 8000, by = 2000)

make_radial <- function(scn) {
  d <- df_long %>% filter(scenario == scn)

  ggplot(d, aes(x = id, y = value, fill = AA)) +
    annotate("segment",
             x = 0.5, xend = 9.5,
             y = breaks_y, yend = breaks_y,
             color = "grey70", linetype = "dotted", linewidth = 0.25) +
    geom_col(width = 0.85, color = "black", linewidth = 0.1) +
    geom_segment(data = rec_df,
                 aes(x = id - 0.42, xend = id + 0.42,
                     y = rec, yend = rec),
                 inherit.aes = FALSE,
                 color = "#333333", linewidth = 0.35) +
    geom_hline(yintercept = ring_y, linewidth = 0.3) +
    scale_fill_manual(values = aa_colors) +
    scale_y_continuous(limits = c(ymin_inner, ymax * 1.2)) +
    scale_x_continuous(limits = c(0.5, 9.5)) +
    coord_polar(start = 0) +
    theme_void() +
    theme(legend.position = "none")
}

size_cm <- 4.5

ggsave(file.path(output_dir, "radial_Open_Flow.png"),
       make_radial("Open_Flow"),
       width = size_cm, height = size_cm,
       units = "cm", bg = "transparent", dpi = 600)

ggsave(file.path(output_dir, "radial_Dense_Array.png"),
       make_radial("Dense_Array"),
       width = size_cm, height = size_cm,
       units = "cm", bg = "transparent", dpi = 600)

ggsave(file.path(output_dir, "radial_Semi_Enclosed.png"),
       make_radial("Semi_Enclosed"),
       width = size_cm, height = size_cm,
       units = "cm", bg = "transparent", dpi = 600)
