# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  Stage 2 — 格网几何计算                                                      ║
# ║                                                                              ║
# ║  目的：由数据实际范围推导一切空间参数，封装为 GridGeometry 结构体。              ║
# ║                                                                              ║
# ║  关键知识：1° 经度 ≠ 固定距离。纬度线平行 → 间距不变；                           ║
# ║  经线向两极汇聚 → 在纬度 φ 处 1° 经度 = 111,320m × cos(φ)。                    ║
# ║                                                                              ║
# ║  resolution 作为参数而非硬编码 → 改 0.5m 格网只需改调用方，函数无需变动。         ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

function ComputeGridGeometry(lons, lats, resolution_m)
    println("[2/6] 计算湖泊边界和格网尺寸...")

    lon_min, lon_max = extrema(lons)
    lat_min, lat_max = extrema(lats)
    mid_lat = (lat_min + lat_max) / 2.0

    # cos 校正：1° 经度在测区纬度处的实际地面距离
    meters_per_deg_lon = METERS_PER_DEG_LAT * cosd(mid_lat)

    # 1m 地面距离 = 多少度
    lon_step = resolution_m / meters_per_deg_lon
    lat_step = resolution_m / METERS_PER_DEG_LAT

    # 测区跨度（米）
    width_m  = (lon_max - lon_min) / lon_step
    height_m = (lat_max - lat_min) / lat_step

    # 向上取整 → 宁可多一列/行也不丢数据
    cols = ceil(Int, width_m)
    rows = ceil(Int, height_m)

    @printf("  经度范围: %.8f → %.8f\n", lon_min, lon_max)
    @printf("  纬度范围: %.8f → %.8f\n", lat_min, lat_max)
    @printf("  中心纬度: %.4f°, cos=%.4f → 1°经度 = %.1f m\n",
            mid_lat, cosd(mid_lat), meters_per_deg_lon)
    @printf("  测区尺寸: %.1f m × %.1f m（%.2f 公顷）\n",
            width_m, height_m, width_m * height_m / 10000)
    @printf("  格网: %d 列 × %d 行 = %d 单元\n", cols, rows, cols * rows)

    return GridGeometry(lon_min, lon_max, lat_min, lat_max, mid_lat,
                        meters_per_deg_lon, lon_step, lat_step,
                        width_m, height_m, cols, rows)
end
