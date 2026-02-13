# ============================================================
# Fig1 (FINAL):
# - Coastline suitability (baseline=blue, cooling-added=red) from SST .nc
# - Mean suitable latitude (dashed) computed from BLUE/RED mask points
# - Yellow site marker sizes & size legend labels are DATA-DRIVEN from SIZE_COL
# - Filters sites by cultivable == True
# - Exports: PNG + PDF + mean-lat table + zoom panels
# - NO area_data; NO override blocks
# ============================================================

import os
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from shapely.geometry import LineString, MultiLineString, GeometryCollection, Polygon, MultiPolygon
from shapely.ops import linemerge
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from cartopy.io import shapereader as shpreader

# ===================== USER SETTINGS =====================

SST_PATH    = "cmems_mod_glo_phy-thetao_anfc_0.083deg_P1M-m_114.nc"
SST_VAR     = "thetao"
LON_NAME    = "longitude"
LAT_NAME    = "latitude"

SPECIES_TEMP_XLS = "species_temp_ranges.xlsx"
SITES_XLS        = "site_species_cover10_needcool_coast20km_with_area.xlsx"
SITES_SHEET      = "Sites_Coast20km"

# columns in SITES_XLS
LON_COL     = "src_lon_deg"
LAT_COL     = "src_lat_deg"
SPECIES_COL = "species"
CULTIVABLE_COL = "cultivable"

# ===== choose ONE area column for marker size =====
# options:
#   "area_Open_Flow_km2"
#   "area_Dense_Array_km2"
#   "area_Semi_Enclosed_km2"
SIZE_COL    = "area_Semi_Enclosed_km2"
AREA_UNIT   = "km²"

# thresholds (months)
BASELINE_THR_MONTHS = 9
ADDED_THR_MONTHS    = 9

DELTA_T = 4.0  # cooling -4C

# nearshore sampling offsets
NEARSHORE_OFFSETS = [
    (0.0, 0.0),
    (0.5, 0.0),
    (-0.5, 0.0),
    (0.0, 0.5),
    (0.0, -0.5),
]

# coastline densify/smooth
MAX_GAP_POINTS   = 8
MIN_SEG_LEN_DEG  = 0.7
DENSIFY_STEP_DEG = 0.15

# output figure
DPI = 600
GRID_FIGSIZE = (11, 16)
FIGSIZE_SINGLE = (10, 12)

USE_FIXED_EXTENT = True
MAP_LON_MIN, MAP_LON_MAX = -130, 150
MAP_LAT_MIN, MAP_LAT_MAX = -80, 80

# colors
COAST_BASELINE_COLOR = "#1f77b4"
COAST_ADDED_COLOR    = "#d62728"
COAST_OTHER_COLOR    = "lightgray"

SITE_OTHER_COLOR     = "0.35"
SITE_FOCAL_FACE      = "#FFD92F"

MEAN_BASELINE_COLOR  = "#9bbf8a"
MEAN_ADDED_COLOR     = "#f79059"

ZOOM_SIZE_SCALE_FACTOR = 100.0

SPECIES_ZOOM_BBOXES = {
    "salmon":     [(-13.0, 10.0, 32.0, 53.0)],
    "cod":        [(-13.0, 10.0, 47.0, 57.0)],
    "grouper":    [(100.0, 113.0, -12.0, -2.0)],
    "barramundi": [(100.0, 113.0, -12.0, -1.0)],
    "seabass":    [(110.0, 125.0, 18.0, 30.0)],
    "seabream":   [(110.0, 135.0, 18.0, 40.0)],
}

# ===================== HELPERS =====================

def densify_line(line, max_step_deg=DENSIFY_STEP_DEG):
    length = line.length
    n = max(int(length / max_step_deg), 1)
    pts = [line.interpolate(float(i) / n, normalized=True) for i in range(n + 1)]
    return pts

def sample_sst_points(lons, lats, sst_da, lon_name=LON_NAME, lat_name=LAT_NAME):
    lons = np.asarray(lons)
    lats = np.asarray(lats)
    data_lons = sst_da[lon_name]

    # handle 0..360 longitude datasets
    if float(data_lons.max()) > 180:
        lons_sample = (lons + 360) % 360
    else:
        lons_sample = lons

    out = sst_da.interp({
        lon_name: ("points", lons_sample),
        lat_name: ("points", lats),
    })
    return out.values

def smooth_mask(mask, max_gap=MAX_GAP_POINTS):
    mask = np.asarray(mask, dtype=bool)
    n = len(mask)
    out = mask.copy()
    i = 0
    while i < n:
        if not mask[i]:
            j = i + 1
            while j < n and not mask[j]:
                j += 1
            gap_len = j - i
            if i > 0 and j < n and gap_len <= max_gap and mask[i-1] and mask[j]:
                out[i:j] = True
            i = j
        else:
            i += 1
    return out

def mask_to_segments(xs, ys, mask):
    xs = np.asarray(xs)
    ys = np.asarray(ys)
    mask = np.asarray(mask)

    segs = []
    n = len(xs)
    i = 0
    while i < n:
        if mask[i]:
            j = i + 1
            while j < n and mask[j]:
                j += 1
            if j - i >= 2:
                coords = [(xs[k], ys[k]) for k in range(i, j)]
                segs.append(LineString(coords))
            i = j
        else:
            i += 1
    return segs

def merge_and_filter_segments(segments, min_len_deg=MIN_SEG_LEN_DEG):
    if not segments:
        return []
    multi = MultiLineString(segments)
    merged_geom = linemerge(multi)

    result = []
    if isinstance(merged_geom, LineString):
        if merged_geom.length >= min_len_deg:
            result.append(merged_geom)
    elif isinstance(merged_geom, MultiLineString):
        for g in merged_geom.geoms:
            if g.length >= min_len_deg:
                result.append(g)
    elif isinstance(merged_geom, GeometryCollection):
        for g in merged_geom.geoms:
            if isinstance(g, (LineString, MultiLineString)):
                if isinstance(g, LineString) and g.length >= min_len_deg:
                    result.append(g)
                elif isinstance(g, MultiLineString):
                    for g2 in g.geoms:
                        if g2.length >= min_len_deg:
                            result.append(g2)
    return result

def plot_lines_from_geom(ax, geom, **kwargs):
    if geom.is_empty:
        return
    if isinstance(geom, LineString):
        x, y = geom.xy
        ax.plot(x, y, **kwargs)
    elif isinstance(geom, MultiLineString):
        for g in geom.geoms:
            x, y = g.xy
            ax.plot(x, y, **kwargs)
    elif isinstance(geom, GeometryCollection):
        for g in geom.geoms:
            if isinstance(g, (LineString, MultiLineString)):
                plot_lines_from_geom(ax, g, **kwargs)

def plot_country_geom(ax, geom, **kwargs):
    if geom is None:
        return
    if isinstance(geom, Polygon):
        x, y = geom.exterior.xy
        ax.plot(x, y, **kwargs)
    elif isinstance(geom, MultiPolygon):
        for poly in geom.geoms:
            x, y = poly.exterior.xy
            ax.plot(x, y, **kwargs)

def coastal_nearshore_masks(xs, ys, theta_da, Tmin, Tmax,
                            delta_t=DELTA_T,
                            lon_name=LON_NAME, lat_name=LAT_NAME,
                            baseline_thr=BASELINE_THR_MONTHS,
                            added_thr=ADDED_THR_MONTHS):
    """
    baseline_ok: months_baseline >= baseline_thr
    added_ok: months_baseline < baseline_thr AND months_union(baseline OR cooled) >= added_thr
    """
    xs = np.asarray(xs)
    ys = np.asarray(ys)
    n  = len(xs)

    if "time" not in theta_da.dims:
        base_any = np.zeros(n, dtype=bool)
        cool_any = np.zeros(n, dtype=bool)
        for dlat, dlon in NEARSHORE_OFFSETS:
            sst_vals = sample_sst_points(xs + dlon, ys + dlat, theta_da, lon_name, lat_name)
            valid = np.isfinite(sst_vals)
            if not valid.any():
                continue
            sst_valid = sst_vals.copy()
            sst_valid[~valid] = np.nan

            base_ok_i = (sst_valid >= Tmin) & (sst_valid <= Tmax)
            cool_ok_i = ((sst_valid - delta_t) >= Tmin) & ((sst_valid - delta_t) <= Tmax)

            base_any |= (base_ok_i & valid)
            cool_any |= (cool_ok_i & valid)

        base_ok  = base_any
        added_ok = cool_any & (~base_ok)
        return base_ok, added_ok

    time_len = theta_da.sizes.get("time")
    base_any_month  = np.zeros((time_len, n), dtype=bool)
    union_any_month = np.zeros((time_len, n), dtype=bool)

    for dlat, dlon in NEARSHORE_OFFSETS:
        sst_vals = sample_sst_points(xs + dlon, ys + dlat, theta_da, lon_name, lat_name)
        valid = np.isfinite(sst_vals)
        if not valid.any():
            continue
        sst_valid = sst_vals.copy()
        sst_valid[~valid] = np.nan

        base_ok_i = (sst_valid >= Tmin) & (sst_valid <= Tmax)
        cool_ok_i = ((sst_valid - delta_t) >= Tmin) & ((sst_valid - delta_t) <= Tmax)
        union_ok_i = base_ok_i | cool_ok_i

        base_any_month  |= (base_ok_i & valid)
        union_any_month |= (union_ok_i & valid)

    base_month_counts  = base_any_month.sum(axis=0)
    union_month_counts = union_any_month.sum(axis=0)

    baseline_ok = base_month_counts >= baseline_thr
    added_ok    = (base_month_counts < baseline_thr) & (union_month_counts >= added_thr)
    return baseline_ok, added_ok

def hemi_mean_lat(lat_list):
    if not lat_list:
        return (None, None)
    arr = np.asarray(lat_list, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return (None, None)
    north = arr[arr >= 0]
    south = arr[arr < 0]
    north_mean = float(north.mean()) if north.size else None
    south_abs_mean = float((-south).mean()) if south.size else None
    return north_mean, south_abs_mean


# ===================== LOAD DATA =====================

df_temp = pd.read_excel(SPECIES_TEMP_XLS)

df_sites_all = pd.read_excel(SITES_XLS, sheet_name=SITES_SHEET)

# cultivable filter (required by you)
if CULTIVABLE_COL in df_sites_all.columns:
    df_sites_all = df_sites_all[df_sites_all[CULTIVABLE_COL] == True].copy()

df_sites_all["_species_lc"] = df_sites_all[SPECIES_COL].astype(str).str.lower()

ds = xr.open_dataset(SST_PATH)
theta = ds[SST_VAR]
theta10 = theta.isel(depth=0) if "depth" in theta.dims else theta
if "depth" in theta10.dims and theta10.sizes.get("depth", 1) == 1:
    theta10 = theta10.isel(depth=0)

# Natural Earth
coast_shp    = shpreader.natural_earth(resolution="110m", category="physical", name="coastline")
land_shp     = shpreader.natural_earth(resolution="110m", category="physical", name="land")
country_shp  = shpreader.natural_earth(resolution="110m", category="cultural", name="admin_0_countries")

coast_geoms   = list(shpreader.Reader(coast_shp).geometries())
land_geoms    = list(shpreader.Reader(land_shp).geometries())
country_geoms = list(shpreader.Reader(country_shp).geometries())

# species list (only those having sites)
species_rows = []
for row in df_temp.itertuples():
    sp_lc = str(row.species).lower()
    df_sp = df_sites_all[df_sites_all["_species_lc"] == sp_lc]
    if df_sp.empty:
        continue
    species_rows.append((row, df_sp))

n_species = len(species_rows)
print(f"Will plot {n_species} species.")


# ===================== SIZE LEGEND (DATA-DRIVEN) =====================

vals_all = pd.to_numeric(df_sites_all.get(SIZE_COL, pd.Series([], dtype=float)), errors="coerce")
vals_all = vals_all[np.isfinite(vals_all)]

size_markers = [50.0, 100.0, 150.0, 200.0]  # scatter s in pt^2

if len(vals_all) > 0:
    q25, q50, q75 = np.quantile(vals_all, [0.25, 0.50, 0.75])

    b1 = vals_all[vals_all <= q25]
    b2 = vals_all[(vals_all > q25) & (vals_all <= q50)]
    b3 = vals_all[(vals_all > q50) & (vals_all <= q75)]
    b4 = vals_all[vals_all > q75]

    reps = [
        float(np.nanmedian(b1)) if len(b1) else np.nan,
        float(np.nanmedian(b2)) if len(b2) else np.nan,
        float(np.nanmedian(b3)) if len(b3) else np.nan,
        float(np.nanmedian(b4)) if len(b4) else np.nan,
    ]
    size_labels = [f"~{v:.1f} {AREA_UNIT}" if np.isfinite(v) else "~NA" for v in reps]
else:
    q25 = q50 = q75 = np.nan
    size_labels = ["~NA"] * 4


# ===================== PLOT GRID =====================

ncols = 2
nrows = int(np.ceil(n_species / ncols))

fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=GRID_FIGSIZE, dpi=DPI)
axes = np.array(axes).reshape(nrows, ncols)
axes_flat = axes.ravel()

legend_handles = [
    mpatches.Patch(color=COAST_BASELINE_COLOR, label="Potential suitability (baseline)"),
    mpatches.Patch(color=COAST_ADDED_COLOR,    label="Incremental potential suitability (with cooling)"),
    Line2D([0], [0], linestyle="--", color=MEAN_BASELINE_COLOR, linewidth=1.5,
           label="Mean suitable latitude (baseline)"),
    Line2D([0], [0], linestyle="--", color=MEAN_ADDED_COLOR, linewidth=1.5,
           label="Mean suitable latitude (with cooling)"),
]
for smark, slabel in zip(size_markers, size_labels):
    legend_handles.append(
        Line2D([0], [0],
               marker="o", color="none",
               markerfacecolor=SITE_FOCAL_FACE,
               markeredgecolor="black",
               markersize=float(np.sqrt(smark)),
               label=slabel)
    )

species_segments = {}
species_points = {}
species_other_points = {}
lat_summary_rows = []

for idx, (row, df_sp) in enumerate(species_rows):
    ax = axes_flat[idx]
    species = row.species
    sp_lc = str(species).lower()
    Tmin = float(row.Tmin_C)
    Tmax = float(row.Tmax_C)

    lon_sites_focal = df_sp[LON_COL].values
    lat_sites_focal = df_sp[LAT_COL].values

    # other sites (remove overlaps by rounded coords)
    df_other = df_sites_all[df_sites_all["_species_lc"] != sp_lc].copy()

    df_sp = df_sp.copy()
    df_sp["_lon_r"] = df_sp[LON_COL].round(3)
    df_sp["_lat_r"] = df_sp[LAT_COL].round(3)
    df_other["_lon_r"] = df_other[LON_COL].round(3)
    df_other["_lat_r"] = df_other[LAT_COL].round(3)

    focal_coords = set(zip(df_sp["_lon_r"], df_sp["_lat_r"]))
    mask_not_overlap = ~df_other[["_lon_r", "_lat_r"]].apply(
        lambda r: (r["_lon_r"], r["_lat_r"]) in focal_coords, axis=1
    )
    df_other = df_other[mask_not_overlap].copy()
    lon_other = df_other[LON_COL].values
    lat_other = df_other[LAT_COL].values

    if USE_FIXED_EXTENT:
        lon_min, lon_max = MAP_LON_MIN, MAP_LON_MAX
        lat_min, lat_max = MAP_LAT_MIN, MAP_LAT_MAX
    else:
        lon_min, lon_max = float(np.nanmin(lon_sites_focal)) - 5, float(np.nanmax(lon_sites_focal)) + 5
        lat_min, lat_max = float(np.nanmin(lat_sites_focal)) - 3, float(np.nanmax(lat_sites_focal)) + 3

    baseline_segments = []
    added_segments = []
    lat_base_pts = []
    lat_add_pts  = []

    for geom in coast_geoms:
        if not isinstance(geom, (LineString, MultiLineString)):
            continue

        lines = [geom] if isinstance(geom, LineString) else list(geom.geoms)

        for line in lines:
            minx, miny, maxx, maxy = line.bounds
            if maxx < lon_min - 0.5 or minx > lon_max + 0.5 or \
               maxy < lat_min - 0.5 or miny > lat_max + 0.5:
                continue

            pts = densify_line(line, max_step_deg=DENSIFY_STEP_DEG)
            xs = [p.x for p in pts]
            ys = [p.y for p in pts]

            base_ok, added_ok = coastal_nearshore_masks(
                xs, ys, theta10, Tmin, Tmax,
                delta_t=DELTA_T,
                lon_name=LON_NAME,
                lat_name=LAT_NAME
            )

            base_ok_sm  = smooth_mask(base_ok,  max_gap=MAX_GAP_POINTS)
            added_ok_sm = smooth_mask(added_ok, max_gap=MAX_GAP_POINTS)

            ys_arr = np.asarray(ys, dtype=float)
            lat_base_pts.extend(ys_arr[base_ok_sm].tolist())
            lat_add_pts.extend(ys_arr[added_ok_sm].tolist())

            baseline_segments.extend(mask_to_segments(xs, ys, base_ok_sm))
            added_segments.extend(mask_to_segments(xs, ys, added_ok_sm))

    baseline_segments = merge_and_filter_segments(baseline_segments, min_len_deg=MIN_SEG_LEN_DEG)
    added_segments    = merge_and_filter_segments(added_segments,    min_len_deg=MIN_SEG_LEN_DEG)

    bn, bs = hemi_mean_lat(lat_base_pts)
    an, as_ = hemi_mean_lat(lat_add_pts)

    lat_summary_rows.append({
        "species": sp_lc,
        "baseline_north_mean_deg": bn,
        "baseline_south_mean_absdeg": bs,
        "cooling_north_mean_deg": an,
        "cooling_south_mean_absdeg": as_,
        "delta_north_deg": (an - bn) if (bn is not None and an is not None) else np.nan,
        "delta_south_absdeg": (as_ - bs) if (bs is not None and as_ is not None) else np.nan,
    })

    ax.set_title(species.capitalize(), fontsize=14, fontweight="bold", loc="left")

    # land fill
    for g in land_geoms:
        try:
            x, y = g.exterior.xy
            ax.fill(x, y, color="#f2f2f2", zorder=0)
        except Exception:
            try:
                for poly in g:
                    x, y = poly.exterior.xy
                    ax.fill(x, y, color="#f2f2f2", zorder=0)
            except Exception:
                pass

    # countries
    for g in country_geoms:
        plot_country_geom(ax, g, color="#c0c0c0", linewidth=0.4, zorder=1.2)

    # background coastline
    for g in coast_geoms:
        plot_lines_from_geom(ax, g, color=COAST_OTHER_COLOR, linewidth=0.3, zorder=1)

    # baseline/added coast masks
    for seg in baseline_segments:
        x, y = seg.xy
        ax.plot(x, y, color=COAST_BASELINE_COLOR, linewidth=1.5, zorder=2,
                solid_joinstyle="round", solid_capstyle="round")
    for seg in added_segments:
        x, y = seg.xy
        ax.plot(x, y, color=COAST_ADDED_COLOR, linewidth=1.8, zorder=3,
                solid_joinstyle="round", solid_capstyle="round")

    # mean-lat dashed
    if bn is not None and np.isfinite(bn):
        ax.axhline(bn, linestyle="--", color=MEAN_BASELINE_COLOR, linewidth=1.2, zorder=2.6)
    if bs is not None and np.isfinite(bs):
        ax.axhline(-bs, linestyle="--", color=MEAN_BASELINE_COLOR, linewidth=1.2, zorder=2.6)
    if an is not None and np.isfinite(an):
        ax.axhline(an, linestyle="--", color=MEAN_ADDED_COLOR, linewidth=1.2, zorder=2.7)
    if as_ is not None and np.isfinite(as_):
        ax.axhline(-as_, linestyle="--", color=MEAN_ADDED_COLOR, linewidth=1.2, zorder=2.7)

    # points
    ax.scatter(lon_other, lat_other, s=10, c=SITE_OTHER_COLOR, alpha=0.5, edgecolors="none", zorder=4)

    # focal sizes (bucketed by global quantiles)
    if SIZE_COL in df_sp.columns:
        vals = pd.to_numeric(df_sp[SIZE_COL], errors="coerce").values
        if np.isfinite(vals).any() and np.isfinite(q25) and np.isfinite(q50) and np.isfinite(q75):
            sizes = np.full(len(vals), 100.0)
            sizes[vals <= q25] = 50.0
            sizes[(vals > q25) & (vals <= q50)] = 100.0
            sizes[(vals > q50) & (vals <= q75)] = 150.0
            sizes[vals > q75] = 200.0
            sizes[~np.isfinite(vals)] = 100.0
        else:
            sizes = np.full(len(df_sp), 100.0)
    else:
        sizes = np.full(len(df_sp), 100.0)

    ax.scatter(lon_sites_focal, lat_sites_focal,
               s=sizes, c=SITE_FOCAL_FACE, edgecolors="black",
               linewidths=0.4, alpha=0.7, zorder=5)

    # zoom boxes
    if sp_lc in SPECIES_ZOOM_BBOXES:
        for (box_lon_min, box_lon_max, box_lat_min, box_lat_max) in SPECIES_ZOOM_BBOXES[sp_lc]:
            rect = mpatches.Rectangle(
                (box_lon_min, box_lat_min),
                box_lon_max - box_lon_min,
                box_lat_max - box_lat_min,
                linewidth=1.5, edgecolor="black", facecolor="none", zorder=10
            )
            ax.add_patch(rect)

    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)
    ax.set_aspect("equal", adjustable="box")
    ax.tick_params(labelsize=10, length=3)

    # caches for zoom
    species_segments[sp_lc] = (baseline_segments, added_segments)
    species_points[sp_lc] = (lon_sites_focal, lat_sites_focal, sizes)
    species_other_points[sp_lc] = (lon_other, lat_other)

# hide unused axes
for j in range(n_species, nrows * ncols):
    axes_flat[j].set_visible(False)

# axis labels layout
for r in range(nrows):
    for c in range(ncols):
        ax = axes[r, c]
        if not ax.get_visible():
            continue
        if c == 0:
            ax.set_ylabel("Latitude", fontsize=12)
        else:
            ax.set_ylabel("")
            ax.set_yticklabels([])
        if r == nrows - 1:
            ax.set_xlabel("Longitude", fontsize=12)
        else:
            ax.set_xlabel("")
            ax.set_xticklabels([])

fig.subplots_adjust(left=0.05, right=0.99, top=0.92, bottom=0.08, wspace=0.04, hspace=0.05)

# legends (3 rows)
fig.legend(handles=legend_handles[:2], loc="lower center", ncol=2, frameon=False,
           fontsize=12, bbox_to_anchor=(0.5, 0.03))
fig.legend(handles=legend_handles[2:4], loc="lower center", ncol=2, frameon=False,
           fontsize=11, bbox_to_anchor=(0.5, 0.01))
fig.legend(handles=legend_handles[4:], loc="lower center", ncol=4, frameon=False,
           fontsize=10, bbox_to_anchor=(0.5, -0.01))

out_png = "Figure1_CoastlineSuitability.png"
out_pdf = "Figure1_CoastlineSuitability.pdf"

fig.savefig(out_png, dpi=DPI, facecolor="white", bbox_inches="tight")
fig.savefig(out_pdf, dpi=DPI, facecolor="white", bbox_inches="tight")
plt.close(fig)
print(f"Saved: {out_png} and {out_pdf}")

# export table: mean lat from masks
df_lat = pd.DataFrame(lat_summary_rows).sort_values("species")
df_lat["delta_north_km"] = df_lat["delta_north_deg"] * 111.0
df_lat["delta_south_km"] = df_lat["delta_south_absdeg"] * 111.0
df_lat.to_excel("species_mean_lat_from_lines.xlsx", index=False)
print("Saved: species_mean_lat_from_lines.xlsx")

# zoom panels
for sp_lc, boxes in SPECIES_ZOOM_BBOXES.items():
    if sp_lc not in species_segments:
        continue

    base_segs, add_segs = species_segments[sp_lc]
    lon_focal, lat_focal, sizes_focal = species_points.get(sp_lc, (np.array([]), np.array([]), np.array([])))
    lon_other_pts, lat_other_pts = species_other_points.get(sp_lc, (np.array([]), np.array([])))

    for idx_box, (box_lon_min, box_lon_max, box_lat_min, box_lat_max) in enumerate(boxes):
        fig_zoom, ax_zoom = plt.subplots(figsize=FIGSIZE_SINGLE, dpi=DPI)

        # thick black border
        for spine in ax_zoom.spines.values():
            spine.set_zorder(9999)
            spine.set_edgecolor("black")
            spine.set_linewidth(20.0)

        # land fill
        for g in land_geoms:
            try:
                x, y = g.exterior.xy
                ax_zoom.fill(x, y, color="#f2f2f2", zorder=0)
            except Exception:
                try:
                    for poly in g:
                        x, y = poly.exterior.xy
                        ax_zoom.fill(x, y, color="#f2f2f2", zorder=0)
                except Exception:
                    pass

        for g in country_geoms:
            plot_country_geom(ax_zoom, g, color="#c0c0c0", linewidth=0.4, zorder=1.2)

        for g in coast_geoms:
            plot_lines_from_geom(ax_zoom, g, color=COAST_OTHER_COLOR, linewidth=0.3, zorder=1)

        for seg in base_segs:
            x, y = seg.xy
            ax_zoom.plot(x, y, color=COAST_BASELINE_COLOR, linewidth=20, zorder=2,
                         solid_joinstyle="round", solid_capstyle="round")
        for seg in add_segs:
            x, y = seg.xy
            ax_zoom.plot(x, y, color=COAST_ADDED_COLOR, linewidth=20, zorder=3,
                         solid_joinstyle="round", solid_capstyle="round")

        # other points
        if lon_other_pts.size > 0:
            mask_other = (lon_other_pts >= box_lon_min) & (lon_other_pts <= box_lon_max) & \
                         (lat_other_pts >= box_lat_min) & (lat_other_pts <= box_lat_max)
            if mask_other.any():
                ax_zoom.scatter(lon_other_pts[mask_other], lat_other_pts[mask_other],
                                s=10, c=SITE_OTHER_COLOR, alpha=0.6, edgecolors="none", zorder=4)

        # focal points
        if lon_focal.size > 0:
            mask_focal = (lon_focal >= box_lon_min) & (lon_focal <= box_lon_max) & \
                         (lat_focal >= box_lat_min) & (lat_focal <= box_lat_max)
            if mask_focal.any():
                sizes_zoom = sizes_focal[mask_focal] * ZOOM_SIZE_SCALE_FACTOR
                ax_zoom.scatter(lon_focal[mask_focal], lat_focal[mask_focal],
                                s=sizes_zoom, c=SITE_FOCAL_FACE,
                                edgecolors="black", linewidths=5, alpha=0.7, zorder=5)

        ax_zoom.set_xlim(box_lon_min, box_lon_max)
        ax_zoom.set_ylim(box_lat_min, box_lat_max)
        ax_zoom.set_aspect("equal", adjustable="box")
        ax_zoom.set_xticks([])
        ax_zoom.set_yticks([])
        ax_zoom.set_title("")
        ax_zoom.set_xlabel("")
        ax_zoom.set_ylabel("")
        fig_zoom.subplots_adjust(left=0, right=1, top=1, bottom=0)

        out_zoom_base = f"Figure1_Zoom_{sp_lc}_{idx_box+1}"
        fig_zoom.savefig(out_zoom_base + ".png", dpi=DPI, facecolor="white",
                         bbox_inches="tight", pad_inches=0)
        fig_zoom.savefig(out_zoom_base + ".pdf", dpi=DPI, facecolor="white",
                         bbox_inches="tight", pad_inches=0)
        plt.close(fig_zoom)
        print(f"Saved: {out_zoom_base}.png/.pdf")
