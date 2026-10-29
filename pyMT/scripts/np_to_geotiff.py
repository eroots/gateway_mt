#!/usr/bin/env python3
"""
Convert a text file of lon/lat/elevation (or lat/lon/elevation) points into
a GeoTIFF raster that can be loaded directly into QGIS.

Usage:
    python np_to_geotiff.py input.xyz output.tif --res 0.001
    python np_to_geotiff.py input.xyz output.tif --res 0.001 --lat-lon-order lat,lon,z

Requirements:
    pip install numpy scipy rasterio --break-system-packages
"""

import argparse
import numpy as np
from scipy.interpolate import griddata
import rasterio
from rasterio.transform import from_origin


def load_points(path, col_order):
    """Load an xyz text file and return lon, lat, z arrays (reordered to lon/lat)."""
    data = np.loadtxt(path)
    if data.shape[1] < 3:
        raise ValueError(f"Expected at least 3 columns, got {data.shape[1]}")

    cols = {name: i for i, name in enumerate(col_order)}
    lon = data[:, cols["lon"]]
    lat = data[:, cols["lat"]]
    z = data[:, cols["z"]]
    return lon, lat, z


def build_grid(lon, lat, z, res, method):
    """Interpolate scattered points onto a regular grid at the given resolution (deg)."""
    lon_min, lon_max = lon.min(), lon.max()
    lat_min, lat_max = lat.min(), lat.max()

    grid_lon = np.arange(lon_min, lon_max + res, res)
    grid_lat = np.arange(lat_min, lat_max + res, res)
    grid_x, grid_y = np.meshgrid(grid_lon, grid_lat)

    grid_z = griddata((lon, lat), z, (grid_x, grid_y), method=method)

    # Flip vertically: raster row 0 = top = max latitude
    grid_z = np.flipud(grid_z)
    top_lat = grid_lat.max()

    return grid_z, lon_min, top_lat


def write_geotiff(path, grid_z, lon_min, top_lat, res, crs, nodata=-9999.0):
    grid_z = np.where(np.isnan(grid_z), nodata, grid_z).astype("float32")
    transform = from_origin(lon_min, top_lat, res, res)

    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=grid_z.shape[0],
        width=grid_z.shape[1],
        count=1,
        dtype=grid_z.dtype,
        crs=crs,
        transform=transform,
        nodata=nodata,
    ) as dst:
        dst.write(grid_z, 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Path to input .xyz text file")
    parser.add_argument("output", help="Path to output .tif raster")
    parser.add_argument("--res", type=float, default=0.001,
                         help="Grid resolution in degrees (default: 0.001, ~100m)")
    parser.add_argument("--method", choices=["linear", "cubic", "nearest"], default="linear",
                         help="Interpolation method (default: linear)")
    parser.add_argument("--col-order", default="lon,lat,z",
                         help="Column order in file as comma list, e.g. 'lat,lon,z' "
                              "if your file is lat/long/elevation (default: lon,lat,z)")
    parser.add_argument("--crs", default="EPSG:4326",
                         help="CRS of the input coordinates (default: EPSG:4326, i.e. WGS84 lat/lon)")
    args = parser.parse_args()

    col_order = [c.strip() for c in args.col_order.split(",")]

    print(f"Loading points from {args.input} ...")
    lon, lat, z = load_points(args.input, col_order)
    print(f"  {len(z)} points, lon range [{lon.min():.5f}, {lon.max():.5f}], "
          f"lat range [{lat.min():.5f}, {lat.max():.5f}]")

    print(f"Gridding at {args.res} deg resolution ({args.method}) ...")
    grid_z, lon_min, top_lat = build_grid(lon, lat, z, args.res, args.method)
    print(f"  Grid shape: {grid_z.shape}")

    print(f"Writing GeoTIFF to {args.output} ...")
    write_geotiff(args.output, grid_z, lon_min, top_lat, args.res, args.crs)

    print("Done. Load the .tif directly into QGIS (Layer > Add Layer > Add Raster Layer).")


if __name__ == "__main__":
    main()