#!/usr/bin/env Rscript
#
# 多波束水深可视化 — 三文件深度图
# 自研方案，不依赖 marmap 包

# ── 入口路径 ──
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg) else "."
ROOT     <- normalizePath(dirname(script_path))
DATA_DIR <- file.path(ROOT, "..", "MBJulia", "result")
OUT_DIR  <- file.path(ROOT, "result")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 水深色阶（20阶，深蓝→浅蓝→绿→黄→橙→红）──
BATHY_COLORS <- colorRampPalette(
  c("#023858", "#045A8D", "#3690C0", "#74A9CF", "#A6BDDB",
    "#D0D1E6", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476",
    "#41AB5D", "#ADDD8E", "#D9F0A3", "#FEE391", "#FEC44F",
    "#FE9929", "#EC7014", "#CC4C02", "#993404", "#D73027")
)(200)

# ── 色标图例 ──
drawColorbar <- function(zrange, title = "Depth (m)") {
  n <- length(BATHY_COLORS)
  image(x = 1, y = seq(0, 1, length.out = n + 1),
        z = matrix(seq(0, 1, length.out = n), nrow = 1),
        col = rev(BATHY_COLORS), axes = FALSE, xlab = "", ylab = "")
  ticks <- pretty(zrange, n = 8)
  tick_pos <- (ticks - zrange[1]) / diff(zrange)
  axis(4, at = tick_pos, labels = ticks, las = 1, cex.axis = 0.7)
  mtext(title, side = 4, line = 2.5, cex = 0.8)
  box()
}

# ── 比例尺 ──
drawScaleMeters <- function(usr, lat_ref, length_m = 50,
                             x = "bottomleft", inset = 5) {
  insetx <- abs(diff(usr[1:2]) * inset / 100)
  insety <- abs(diff(usr[3:4]) * inset / 100)
  X <- switch(x,
    bottomright = usr[2] - insetx, topright = usr[2] - insetx,
    bottomleft  = usr[1] + insetx, topleft  = usr[1] + insetx)
  Y <- switch(x,
    bottomright = usr[3] + insety, topright = usr[4] - insety,
    bottomleft  = usr[3] + insety, topleft  = usr[4] - insety)
  cos_lat <- cos(2 * pi * lat_ref / 360)
  perdeg <- 2 * pi * (6372.798 + 21.38 * cos_lat) * cos_lat / 360
  deg_len <- length_m / (perdeg * 1000)
  arrows(X, Y, X + deg_len, Y, code = 3, length = 0.05, angle = 90)
  text(X + deg_len/2, Y, labels = paste(length_m, "m"),
       adj = c(0.5, -0.5), cex = 0.7)
}

# ── 构建格网矩阵 ──
BuildGridMatrix <- function(lons, lats, depths) {
  lon_min <- min(lons); lon_max <- max(lons)
  lat_min <- min(lats); lat_max <- max(lats)
  mid_lat <- (lat_min + lat_max) / 2
  m_per_deg_lon <- 111320.0 * cos(mid_lat * pi / 180)
  lon_step <- 1.0 / m_per_deg_lon
  lat_step <- 1.0 / 111320.0
  cols <- round((lon_max - lon_min) / lon_step) + 1
  rows <- round((lat_max - lat_min) / lat_step) + 1

  mat <- matrix(NA_real_, nrow = rows, ncol = cols)
  for (i in seq_along(lons)) {
    c <- round((lons[i] - lon_min) / lon_step) + 1
    r <- round((lats[i] - lat_min) / lat_step) + 1
    if (c >= 1 && c <= cols && r >= 1 && r <= rows) {
      mat[r, c] <- depths[i]
    }
  }
  lon_edges <- seq(lon_min - lon_step/2, by = lon_step, length.out = cols + 1)
  lat_edges <- seq(lat_min - lat_step/2, by = lat_step, length.out = rows + 1)

  list(mat = mat, lon_edges = lon_edges, lat_edges = lat_edges, mid_lat = mid_lat)
}

# ── 单幅深度图 ──
DepthPlot <- function(data_file, out_png, title) {
  d <- read.table(data_file, col.names = c("lon", "lat", "depth"))
  cat(sprintf("  %s: %d 点, %.2f ~ %.2f m\n",
              basename(data_file), nrow(d), min(d$depth), max(d$depth)))

  g <- BuildGridMatrix(d$lon, d$lat, d$depth)
  zrange <- range(g$mat, na.rm = TRUE)
  dr <- diff(zrange)
  step_major <- if (dr > 20) 2 else if (dr > 10) 1 else 0.5

  png(out_png, width = 14, height = 12, units = "in", res = 250)
  layout(matrix(c(1, 2), 1, 2), widths = c(7, 1))
  par(mar = c(4, 4, 3, 0.5))

  image(g$lon_edges, g$lat_edges, t(g$mat),
        col = rev(BATHY_COLORS), zlim = zrange,
        xlab = "Longitude (\u00B0E)", ylab = "Latitude (\u00B0N)",
        main = "", useRaster = TRUE)
  box()

  # 主等深线（粗线 + 标注）
  contour(g$lon_edges[-1] - diff(g$lon_edges)[1]/2,
          g$lat_edges[-1] - diff(g$lat_edges)[1]/2,
          t(g$mat),
          levels = seq(floor(zrange[1]), ceiling(zrange[2]),
                       by = step_major),
          lwd = 0.6, lty = 1, col = "grey20",
          drawlabels = TRUE, add = TRUE)

  # 辅等深线（细虚线，半间距）
  contour(g$lon_edges[-1] - diff(g$lon_edges)[1]/2,
          g$lat_edges[-1] - diff(g$lat_edges)[1]/2,
          t(g$mat),
          levels = seq(floor(zrange[1]), ceiling(zrange[2]),
                       by = step_major / 2),
          lwd = 0.2, lty = 2, col = "grey50",
          drawlabels = FALSE, add = TRUE)

  drawScaleMeters(par("usr"), g$mid_lat, length_m = 50)

  par(mar = c(4, 1, 3, 3))
  drawColorbar(zrange)
  dev.off()
  cat(sprintf("    \u2192 %s\n", out_png))
}

# ═══════════════════════════════════════════════
cat("========================================\n")
cat("多波束水深可视化\n")
cat("========================================\n\n")

targets <- list(
  list(f = "gridded_1m_xyz.txt",
       t = "Lake Bathymetry \u2014 MBSystem Grid Average"),
  list(f = "gridded_1m_clean.txt",
       t = "Lake Bathymetry \u2014 Julia Filtered"),
  list(f = "gridded_1m_interpolated.txt",
       t = "Lake Bathymetry \u2014 Julia Interpolated")
)

for (x in targets) {
  fp <- file.path(DATA_DIR, x$f)
  if (file.exists(fp)) {
    DepthPlot(fp, file.path(OUT_DIR, sub("\\.txt$", ".png", x$f)), x$t)
  } else {
    cat(sprintf("  \u26A0 跳过 %s（文件不存在）\n", x$f))
  }
}

cat("\n完成\n")
