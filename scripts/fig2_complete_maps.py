import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Wedge, Circle
from matplotlib.lines import Line2D
from shapely.ops import unary_union
import geopandas as gpd


# ==========================================
# 1. Global settings
# ==========================================
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.size"] = 8
plt.rcParams["figure.dpi"] = 300


# ==========================================
# 2. Colors
# ==========================================
aa_colors = {
    "Ile": "#41b6c4", "Leu": "#225ea8", "Val": "#a1dab4",
    "AAA": "#fecc5c", "Trp": "#fd8d3c",
    "His": "#9e9ac8", "Lys": "#54278f",
    "Thr": "#f03b20", "SAA": "#bdbdbd"
}
aa_order = ["Ile", "Leu", "Val", "AAA", "Trp", "His", "Lys", "Thr", "SAA"]

pie_color_high = "#FFD92F"  # SSR ≥ 90%
pie_color_low = "#8E4AB9"   # SSR < 90%


# ==========================================
# 3. Helpers
# ==========================================
def _coerce_ratio(x):
    """Convert SSR to [0,1]. Supports 0-1 or 0-100 (%) scales."""
    if pd.isna(x):
        return np.nan
    try:
        v = float(x)
    except Exception:
        return np.nan
    if v > 1.5:
        v = v / 100.0
    return v


def load_ssr_map(ssr_file, scenario_key):
    """
    Load SSR mapping (country -> ssr) from any sheet.
    Assumes columns are named:
      SSR_Open_Flow, SSR_Dense_Array, SSR_Semi_Enclosed
    """
    target_col = f"SSR_{scenario_key}"

    xls = pd.ExcelFile(ssr_file)
    for sh in xls.sheet_names:
        df = pd.read_excel(ssr_file, sheet_name=sh)

        if "Country" not in df.columns and "country" in df.columns:
            df = df.rename(columns={"country": "Country"})
        if "Country" not in df.columns:
            continue
        if target_col not in df.columns:
            continue

        out = df[["Country", target_col]].copy()
        out[target_col] = out[target_col].apply(_coerce_ratio)
        return dict(zip(out["Country"], out[target_col]))

    raise ValueError(f"Could not find column '{target_col}' in any sheet of {ssr_file}.")


def prepare_world(shapefile_path):
    world = gpd.read_file(shapefile_path)

    col = "ADMIN" if "ADMIN" in world.columns else ("NAME_EN" if "NAME_EN" in world.columns else "name")
    world = world.rename(columns={col: "name"})
    world["original_name"] = world["name"]

    name_map = {
        "United States of America": "United States",
        "Dem. Rep. Congo": "Congo",
        "Dominican Rep.": "Dominican Republic",
        "South Korea": "Korea, Republic of",
        "North Korea": "Korea, Democratic People's Republic of",
        "Czechia": "Czech Republic",
        "Myanmar": "Burma",
        "Taiwan": "China",
    }
    world["name"] = world["name"].replace(name_map)

    # Merge Taiwan geometry into China if present
    if "China" in world["name"].values and "Taiwan" in world["original_name"].values:
        china_geom = unary_union(world.loc[world["name"].isin(["China", "Taiwan"]), "geometry"])
        world = world[~world["name"].isin(["China", "Taiwan"])]
        world = pd.concat(
            [world, gpd.GeoDataFrame([{"name": "China", "geometry": china_geom}], crs=world.crs)],
            ignore_index=True,
        )

    return world[world["name"] != "Antarctica"]


def read_coverage_df(excel_path, scenario_key):
    """
    Read coverage sheet named:
      coverage_Open_Flow, coverage_Dense_Array, coverage_Semi_Enclosed
    and expects a 'TotalEAA' column (case-insensitive).
    """
    sheet_name = f"coverage_{scenario_key}"
    df = pd.read_excel(excel_path, sheet_name=sheet_name)

    if "Country" in df.columns:
        df = df.rename(columns={"Country": "country"})
    elif "country" not in df.columns:
        raise ValueError(f"Missing Country/country column in sheet: {sheet_name}")

    # TotalEAA (case-insensitive)
    val_col = next((c for c in df.columns if str(c).strip().lower() == "totaleaa"), None)
    if val_col is None:
        raise ValueError(f"Missing TotalEAA column in sheet: {sheet_name}")

    df["coverage"] = pd.to_numeric(df[val_col], errors="coerce")
    return df[["country", "coverage"]]


# ==========================================
# 4. Plot panel
# ==========================================
def plot_panel(
    ax,
    coverage_df,
    title,
    panel_label,
    ssr_map,
    world,
    cmap,
    norm,
    lon_range,
    lat_range,
    ssr_color_threshold=0.90,
    min_pie_radius_factor=0.015,
    max_pie_radius_factor=0.035,
):
    merged = world.merge(coverage_df, how="left", left_on="name", right_on="country")
    merged["value"] = merged["coverage"].clip(lower=norm.vmin, upper=norm.vmax)

    ax.axis("off")
    merged.boundary.plot(ax=ax, linewidth=0.2, color="#404040")
    merged.plot(column="value", ax=ax, cmap=cmap, norm=norm)

    ax.set_xlim(lon_range)
    ax.set_ylim(lat_range)

    ax.set_title(title, fontsize=12, fontweight="bold", pad=6)
    ax.text(
        0.0,
        1.0,
        panel_label,
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        va="bottom",
        ha="right",
    )

    # Scale pie size by bbox area proxy
    base_range = min(lon_range[1] - lon_range[0], lat_range[1] - lat_range[0])
    geoms = merged["geometry"]

    areas = geoms.apply(
        lambda g: (g.bounds[2] - g.bounds[0]) * (g.bounds[3] - g.bounds[1])
        if g is not None and (not g.is_empty) else 0
    )
    sqrt_areas = np.sqrt(areas)
    valid = sqrt_areas[sqrt_areas > 0]
    if len(valid) == 0:
        return

    s_min, s_max = valid.min(), valid.max()
    denom = (s_max - s_min) if (s_max != s_min) else 1.0
    min_r = base_range * min_pie_radius_factor
    max_r = base_range * max_pie_radius_factor

    for idx, row in merged.iterrows():
        # Taiwan merged into China
        if row.get("original_name") == "Taiwan":
            continue

        name = row["name"]
        if name not in ssr_map:
            continue

        ssr = _coerce_ratio(ssr_map[name])
        if pd.isna(ssr):
            continue

        geom = row["geometry"]
        if geom is None or geom.is_empty:
            continue

        s_area = sqrt_areas.iloc[idx]
        rad = min_r + (s_area - s_min) / denom * (max_r - min_r)
        rad = float(np.clip(rad, min_r, max_r))

        cx, cy = geom.centroid.x, geom.centroid.y
        angle = 360.0 * float(np.clip(ssr, 0, 1))
        color = pie_color_high if ssr >= ssr_color_threshold else pie_color_low

        ax.add_patch(Wedge((cx, cy), rad, 90, 90 + angle, facecolor=color, edgecolor=None, zorder=10))
        ax.add_patch(Wedge((cx, cy), rad, 90 + angle, 450, facecolor="white", edgecolor=None, zorder=10))
        ax.add_patch(Circle((cx, cy), rad, facecolor="none", edgecolor="black", linewidth=0.3, zorder=11))


# ==========================================
# 5. Build full figure
# ==========================================
def make_complete_pdf(
    excel_path,
    shapefile_path,
    ssr_file,
    output_path,
    lon_range=(-170, 160),
    lat_range=(-55, 85),
):
    world = prepare_world(shapefile_path)

    cmap = mcolors.LinearSegmentedColormap.from_list(
        "cov",
        [
            (0.0, "#dbeaf7"),
            (0.25, "#87b6d9"),
            (0.499, "#1f78b4"),
            (0.501, "#f4b5b5"),
            (0.75, "#e77b71"),
            (1.0, "#a62019"),
        ],
        N=256,
    )

    vmin, vmax = 0.0, 2.0
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=1.0, vmax=vmax)

    # ---- strict scenario naming (no compatibility) ----
    scenarios = ["Open_Flow", "Dense_Array", "Semi_Enclosed"]

    df_open = read_coverage_df(excel_path, "Open_Flow")
    df_dense = read_coverage_df(excel_path, "Dense_Array")
    df_semi = read_coverage_df(excel_path, "Semi_Enclosed")

    ssr_open = load_ssr_map(ssr_file, "Open_Flow")
    ssr_dense = load_ssr_map(ssr_file, "Dense_Array")
    ssr_semi = load_ssr_map(ssr_file, "Semi_Enclosed")

    fig = plt.figure(figsize=(8, 12))
    gs = fig.add_gridspec(3, 1, hspace=0.1, left=0.05, right=0.95, top=0.98, bottom=0.18)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])

    plot_panel(ax1, df_open, "Open_Flow", "a", ssr_open, world, cmap, norm, lon_range, lat_range)
    plot_panel(ax2, df_dense, "Dense_Array", "b", ssr_dense, world, cmap, norm, lon_range, lat_range)
    plot_panel(ax3, df_semi, "Semi_Enclosed", "c", ssr_semi, world, cmap, norm, lon_range, lat_range)

    # --------------------------
    # Bottom legends
    # --------------------------
    # 1) Colorbar
    cax = fig.add_axes([0.25, 0.12, 0.5, 0.012])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")

    ticks = [0.5, 1.0, 1.5, 2.0]
    tick_labels = [f">{int(vmax * 100)}%" if t == vmax else f"{int(t * 100)}%" for t in ticks]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(tick_labels)
    cbar.set_label("Coverage (supply / recommended)", fontsize=11, labelpad=5)
    cbar.outline.set_linewidth(0.5)
    cbar.ax.tick_params(labelsize=10, width=0.5)

    # 2) SSR legend (circles)
    ssr_legend_elements = [
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor=pie_color_high, markeredgecolor="black",
               markeredgewidth=0.5, markersize=8, label="SSR ≥ 90%"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor=pie_color_low, markeredgecolor="black",
               markeredgewidth=0.5, markersize=8, label="SSR < 90%"),
    ]
    fig.legend(
        handles=ssr_legend_elements,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.08),
        ncol=2,
        frameon=False,
        fontsize=10,
    )

    # 3) Amino acid legend (squares)
    aa_legend_elements = [
        Line2D([0], [0], marker="s", color="none",
               markerfacecolor=aa_colors[aa],
               markeredgecolor="white", markeredgewidth=0.5,
               markersize=8, label=aa)
        for aa in aa_order
    ]
    fig.legend(
        handles=aa_legend_elements,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.04),
        ncol=9,
        frameon=False,
        fontsize=10,
        handletextpad=0.2,
        columnspacing=0.8,
    )

    out_dir = os.path.dirname(output_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    fig.savefig(output_path, bbox_inches="tight", format="pdf")
    plt.close(fig)


if __name__ == "__main__":
    # Repo-friendly paths
    excel_file = "data/input/country_eaa_weighted_recommendations_and_coverage_with_total.xlsx"
    ssr_file = "data/input/SSR_IDR_SSRcapped_output.xlsx"
    shapefile = "data/shapefile/ne_50m_admin_0_countries.shp"

    output_file = "results/Figure2_Complete_Squares.pdf"

    make_complete_pdf(excel_file, shapefile, ssr_file, output_file)
