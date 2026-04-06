#!/usr/bin/env python3
import numpy as np
import argparse

def load_data(file_path, skip_header=True):
    return np.loadtxt(file_path, skiprows=1 if skip_header else 0, dtype=str)

def filter_by_radius(data, x_index, y_index, x0, y0, radius):
    x = data[:, x_index].astype(float)
    y = data[:, y_index].astype(float)
    distances = np.sqrt((x - x0)**2 + (y - y0)**2)
    return data[distances <= radius]

def save_filtered_data(output_path, filtered_data):
    np.savetxt(output_path, filtered_data, fmt='%s')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter rows by distance from (x0, y0)")
    parser.add_argument("--input", required=True, help="Input data file")
    parser.add_argument("--output", required=True, help="Output filtered data file")
    parser.add_argument("--x0", type=float, required=True, help="X center coordinate")
    parser.add_argument("--y0", type=float, required=True, help="Y center coordinate")
    parser.add_argument("--radius", type=float, required=True, help="Radius to include")
    parser.add_argument("--xcol", type=int, required=True, help="Index of x column (0-based)")
    parser.add_argument("--ycol", type=int, required=True, help="Index of y column (0-based)")
    parser.add_argument("--header", action="store_true", help="Set this if the file has a header row")
    args = parser.parse_args()

    data = load_data(args.input, skip_header=args.header)
    filtered = filter_by_radius(data, args.xcol, args.ycol, args.x0, args.y0, args.radius)
    save_filtered_data(args.output, filtered)

    print(f"Filtered {len(filtered)} rows within radius {args.radius} of ({args.x0}, {args.y0})")

'''
Example use (in terminal):
python cut_ogle_map.py \
  --input blg157.5.map \
  --output filtered_blg157.5.map \
  --x0 626.11 --y0 295.21 --radius 650 \
  --xcol 3 --ycol 4 \
  --header
  '''