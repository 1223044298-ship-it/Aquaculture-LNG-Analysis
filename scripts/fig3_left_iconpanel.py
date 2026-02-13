import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import matplotlib.colors as mcolors

# ==========================================
# Global settings
# ==========================================
ICON_OUTLINE_THICKNESS = 4
plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.linewidth"] = 0

# ==========================================
# Utility functions
# ==========================================
def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip("#")
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))


def create_tinted_icon_array(icon_arr, color_segments, outline_thickness=3):
    base_arr = np.array(icon_arr)
    h, w, _ = base_arr.shape
    tinted_arr = np.zeros_like(base_arr)
    tinted_arr[:, :, 3] = base_arr[:, :, 3]

    x_start = 0
    for rgb_color, frac in color_segments:
        if frac <= 0:
            continue
        x_end = x_start + int(round(frac * w))
        x_end = min(x_end, w)
        tinted_arr[:, x_start:x_end, :3] = rgb_color
        x_start = x_end

    if x_start < w and color_segments:
        tinted_arr[:, x_start:, :3] = color_segments[-1][0]

    # outline
    alpha_mask = base_arr[:, :, 3] > 0
    eroded_mask = alpha_mask.copy()
    for dy in [-1, 0, 1]:
        for dx in [-1, 0, 1]:
            if dx == 0 and dy == 0:
                continue
            shifted = np.roll(alpha_mask, shift=(dy, dx), axis=(0, 1))
            eroded_mask &= shifted

    border_mask = alpha_mask & (~eroded_mask)

    if outline_thickness > 1:
        dilated_mask = border_mask.copy()
        neighbor_shifts = [(-1, -1), (-1, 0), (-1, 1),
                           (0, -1), (0, 1),
                           (1, -1), (1, 0), (1, 1)]
        for _ in range(outline_thickness - 1):
            new_mask = dilated_mask.copy()
            for dy, dx in neighbor_shifts:
                shifted = np.roll(border_mask, shift=(dy, dx), axis=(0, 1))
                new_mask |= shifted
            dilated_mask = new_mask
        border_mask = dilated_mask

    tinted_arr[border_mask, :3] = (0, 0, 0)
    tinted_arr[border_mask, 3] = base_arr[border_mask, 3]

    return tinted_arr


def build_icon_segments(values, colors_rgb, units_per_icon=50):
    icons = []
    current_icon = []
    remaining_capacity = 1.0

    units_list = [val / units_per_icon for val in values]

    for value_units, color_rgb in zip(units_list, colors_rgb):
        remaining_units = value_units
        while remaining_units > 1e-9:
            fill_amount = min(remaining_units, remaining_capacity)
            if fill_amount > 1e-9:
                current_icon.append((color_rgb, fill_amount))
                remaining_units -= fill_amount
                remaining_capacity -= fill_amount
            if remaining_capacity <= 1e-9:
                icons.append(current_icon)
                current_icon = []
                remaining_capacity = 1.0

    if current_icon:
        remaining_frac = max(0.0, 1.0 - sum(frac for _, frac in current_icon))
        if remaining_frac > 1e-9:
            current_icon.append(((255, 255, 255), remaining_frac))
        icons.append(current_icon)

    return icons


# ==========================================
# Main plotting function
# ==========================================
def plot_icon_panel(df, icon_path, colors_hex, output_path,
                    top_n=15, units_per_icon=50):

    # required columns
    required_cols = ["Country", "through_k", "flux_k", "weak_k"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")

    df = df.copy()

    for col in ["through_k", "flux_k", "weak_k"]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)

    df["seg_through"] = df["through_k"].clip(lower=0)
    df["seg_flux_add"] = (df["flux_k"] - df["through_k"]).clip(lower=0)
    df["seg_weak_add"] = (df["weak_k"] - df["flux_k"]).clip(lower=0)
    df["total"] = df[["seg_through", "seg_flux_add", "seg_weak_add"]].sum(axis=1)

    df = df.sort_values("total", ascending=False).reset_index(drop=True)

    top_df = df.head(top_n)
    rest_df = df.iloc[top_n:]

    if len(rest_df) > 0:
        others = pd.DataFrame({
            "Country": ["Others"],
            "seg_through": [rest_df["seg_through"].sum()],
            "seg_flux_add": [rest_df["seg_flux_add"].sum()],
            "seg_weak_add": [rest_df["seg_weak_add"].sum()],
            "total": [rest_df["total"].sum()],
        })
        plot_df = pd.concat([top_df, others], ignore_index=True)
    else:
        plot_df = top_df.copy()

    icon_img = Image.open(icon_path).convert("RGBA")
    icon_img = icon_img.rotate(180, expand=True)
    icon_arr = np.array(icon_img)
    colors_rgb = [hex_to_rgb(c) for c in colors_hex]

    icons_per_row = []
    max_icons_count = 0

    for _, row in plot_df.iterrows():
        segments = [row["seg_through"], row["seg_flux_add"], row["seg_weak_add"]]
        icon_list = build_icon_segments(segments, colors_rgb, units_per_icon)
        icons_per_row.append(icon_list)
        max_icons_count = max(max_icons_count, len(icon_list))

    num_rows = len(plot_df)

    fig_height = num_rows * 0.7 + 1.0
    fig_width = (max_icons_count * 0.85) + 1.0

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    max_x_end = 0

    for i, icon_list in enumerate(icons_per_row):
        x_offset = 0
        for icon_segments in icon_list:
            tinted_icon = create_tinted_icon_array(
                icon_arr, icon_segments,
                outline_thickness=ICON_OUTLINE_THICKNESS
            )
            ax.imshow(
                tinted_icon,
                extent=(x_offset, x_offset + 0.8, i - 0.35, i + 0.35),
                zorder=10
            )
            x_offset += 0.85

        max_x_end = max(max_x_end, x_offset)

    ax.set_xlim(0, max_x_end + 0.1)
    ax.set_ylim(num_rows - 0.5, -0.5)
    ax.axis("off")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    fig.savefig(output_path + ".pdf", dpi=600, bbox_inches="tight", pad_inches=0)
    fig.savefig(output_path + ".png", dpi=600, bbox_inches="tight", pad_inches=0, transparent=True)
    plt.close(fig)


# ==========================================
# Entry point
# ==========================================
if __name__ == "__main__":

    input_file = "data/input/child.xlsx"
    icon_file = "assets/child_icon_centered.png"
    output_file = "results/Figure3_Left_IconPanel"

    df = pd.read_excel(input_file)

    colors_hex = ["#2A9D8F", "#E9C46A", "#E76F51"]

    plot_icon_panel(df, icon_file, colors_hex, output_file)
