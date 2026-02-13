import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors

# ---------------------------
# Global config
# ---------------------------
MAIN_FIG_SIZE_CM = 8.5
MAIN_FIG_SIZE_INCH = MAIN_FIG_SIZE_CM / 2.54

LINE_WIDTH = 2.0
FONT_FAMILY = "Arial"

plt.rcParams["font.family"] = FONT_FAMILY
plt.rcParams["axes.linewidth"] = 0

# ---------------------------
# Colors & styles
# ---------------------------
stroke_colors = {
    "Fish": "#2E7D32",
    "Beef": "#1565C0",
}

line_styles = {
    "Fish": "--",  # dashed
    "Beef": "-",   # solid
}

new_red_colors = [
    "#FFF176", "#FFEE58", "#FFCA28", "#FFB300", "#FF9800",
    "#F57C00", "#FF5722", "#E64A19", "#D32F2F", "#B71C1C"
]
cmap_fill = mcolors.LinearSegmentedColormap.from_list("single_hue_red", new_red_colors)

USD_LIMIT_MIN = 150
USD_LIMIT_MAX = 750
CO2_LIMIT_MAX = 5000
POWER_FACTOR = 0.55
BASE_SCALE_FACTOR = 0.03
COLOR_POWER_FACTOR = 1.0

def get_fill_color(value):
    clamped_val = min(max(value, USD_LIMIT_MIN), USD_LIMIT_MAX)
    norm = (clamped_val - USD_LIMIT_MIN) / (USD_LIMIT_MAX - USD_LIMIT_MIN)
    return cmap_fill(norm ** COLOR_POWER_FACTOR)

def calculate_radius(val):
    clamped_val = min(val, CO2_LIMIT_MAX)
    return np.sqrt(clamped_val ** POWER_FACTOR) * BASE_SCALE_FACTOR

# ---------------------------
# Data (static for reproducibility)
# ---------------------------
data_co2 = {
    "Lys": {"Fish": 179.0,  "Beef": 1530.6},
    "SAA": {"Fish": 406.5,  "Beef": 3933.9},
    "His": {"Fish": 574.7,  "Beef": 3908.3},
    "Ile": {"Fish": 355.5,  "Beef": 3149.3},
    "Leu": {"Fish": 202.3,  "Beef": 1700.9},
    "AAA": {"Fish": 224.6,  "Beef": 1895.6},
    "Thr": {"Fish": 369.7,  "Beef": 3108.8},
    "Trp": {"Fish": 1480.0, "Beef": 13111.9},
    "Val": {"Fish": 321.4,  "Beef": 2930.3},
}

data_usd = {
    "Lys": {"Fish": 214.4,  "Beef": 255.1},
    "SAA": {"Fish": 486.7,  "Beef": 655.7},
    "His": {"Fish": 688.3,  "Beef": 651.4},
    "Ile": {"Fish": 425.7,  "Beef": 524.9},
    "Leu": {"Fish": 242.2,  "Beef": 283.5},
    "AAA": {"Fish": 269.0,  "Beef": 315.9},
    "Thr": {"Fish": 442.7,  "Beef": 518.1},
    "Trp": {"Fish": 1772.4, "Beef": 2185.3},
    "Val": {"Fish": 384.8,  "Beef": 488.4},
}

amino_acids = ["Lys", "SAA", "His", "Ile", "Leu", "AAA", "Thr", "Trp", "Val"]

# ---------------------------
# Output directory (repo-relative)
# ---------------------------
OUT_DIR = os.path.join("results", "fig3_right")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------
# Main grid
# ---------------------------
def plot_grid(output_stem: str):
    out_pdf = os.path.join(OUT_DIR, output_stem + ".pdf")
    out_png = os.path.join(OUT_DIR, output_stem + ".png")

    fig, axes = plt.subplots(3, 3, figsize=(MAIN_FIG_SIZE_INCH, MAIN_FIG_SIZE_INCH))
    plt.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99, wspace=0.0, hspace=0.0)
    axes_flat = axes.flatten()

    for i, ax in enumerate(axes_flat):
        if i >= len(amino_acids):
            ax.axis("off")
            continue

        amino = amino_acids[i]
        items = []
        for animal in ["Fish", "Beef"]:
            co2 = data_co2[amino][animal]
            usd = data_usd[amino][animal]
            radius = calculate_radius(co2)
            fill_c = get_fill_color(usd)
            edge_c = stroke_colors[animal]
            l_style = line_styles[animal]
            items.append({"r": radius, "fill": fill_c, "edge": edge_c, "ls": l_style})

        items.sort(key=lambda x: x["r"], reverse=True)

        x_center = 0.5
        base_y = 0.05
        for item in items:
            circle = mpatches.Circle(
                (x_center, base_y + item["r"]),
                item["r"],
                facecolor=item["fill"],
                edgecolor=item["edge"],
                linestyle=item["ls"],
                linewidth=LINE_WIDTH,
                zorder=10,
            )
            ax.add_patch(circle)

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")

    fig.savefig(out_pdf, dpi=600)
    fig.savefig(out_png, dpi=600, transparent=True)
    plt.close(fig)

# ---------------------------
# Legends
# ---------------------------
def plot_legend_species(output_stem: str):
    out_pdf = os.path.join(OUT_DIR, output_stem + ".pdf")
    out_png = os.path.join(OUT_DIR, output_stem + ".png")

    W_inch, H_inch = 1.5, 3.0
    inch_per_unit = MAIN_FIG_SIZE_INCH / 3.0
    x_limit = W_inch / inch_per_unit
    y_limit = H_inch / inch_per_unit

    fig, ax = plt.subplots(figsize=(W_inch, H_inch))
    ax.set_xlim(0, x_limit)
    ax.set_ylim(0, y_limit)
    ax.axis("off")

    r_common = calculate_radius(2500)
    x_center = x_limit * 0.5
    y_center = y_limit * 0.5
    gap = 0.2

    c_beef = mpatches.Circle(
        (x_center, y_center + r_common + gap / 2),
        r_common,
        facecolor="none",
        edgecolor=stroke_colors["Beef"],
        linestyle="-",
        linewidth=LINE_WIDTH,
    )
    ax.add_patch(c_beef)

    c_fish = mpatches.Circle(
        (x_center, y_center - r_common - gap / 2),
        r_common,
        facecolor="none",
        edgecolor=stroke_colors["Fish"],
        linestyle="--",
        linewidth=LINE_WIDTH,
    )
    ax.add_patch(c_fish)

    fig.savefig(out_pdf, dpi=600, bbox_inches="tight")
    fig.savefig(out_png, dpi=600, bbox_inches="tight", transparent=True)
    plt.close(fig)

def plot_legend_carbon(output_stem: str):
    out_pdf = os.path.join(OUT_DIR, output_stem + ".pdf")
    out_png = os.path.join(OUT_DIR, output_stem + ".png")

    W_inch, H_inch = 2.0, 4.0
    inch_per_unit = MAIN_FIG_SIZE_INCH / 3.0
    x_limit = W_inch / inch_per_unit
    y_limit = H_inch / inch_per_unit

    fig, ax = plt.subplots(figsize=(W_inch, H_inch))
    ax.set_xlim(0, x_limit)
    ax.set_ylim(0, y_limit)
    ax.axis("off")

    sizes = [200, 1500, 3000, 5000]
    x_center = x_limit * 0.5
    current_y = y_limit * 0.9

    for s in sizes:
        r = calculate_radius(s)
        c = mpatches.Circle((x_center, current_y - r), r, facecolor="black", edgecolor="none")
        ax.add_patch(c)
        current_y -= (2 * r + 0.15)

    fig.savefig(out_pdf, dpi=600, bbox_inches="tight")
    fig.savefig(out_png, dpi=600, bbox_inches="tight", transparent=True)
    plt.close(fig)

def plot_legend_cost(output_stem: str):
    out_pdf = os.path.join(OUT_DIR, output_stem + ".pdf")
    out_png = os.path.join(OUT_DIR, output_stem + ".png")

    W_inch, H_inch = 3.0, 1.0
    fig, ax = plt.subplots(figsize=(W_inch, H_inch))
    ax.axis("off")

    ax_cbar = ax.inset_axes([0.1, 0.4, 0.8, 0.3])
    gradient = np.linspace(0, 1, 256).reshape(1, -1)
    ax_cbar.imshow(gradient, aspect="auto", cmap=cmap_fill)
    ax_cbar.set_axis_off()

    fig.savefig(out_pdf, dpi=600, bbox_inches="tight")
    fig.savefig(out_png, dpi=600, bbox_inches="tight", transparent=True)
    plt.close(fig)

if __name__ == "__main__":
    plot_grid("Figure3_RightPanel_Main")
    plot_legend_species("Figure3_RightPanel_Legend_Species")
    plot_legend_carbon("Figure3_RightPanel_Legend_Carbon")
    plot_legend_cost("Figure3_RightPanel_Legend_Cost")
