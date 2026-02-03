
# Usage
# python vel2fits.py RT_vel_Nrotate25_+z_HCO+10.dat -ipix=256 -channels=101 | (for one file)
# python vel2fits.py "RT_vel_Nrotate*_+z_*.dat" -ipix=256 -channels=101     | (for multiple files)
# python vel2fits.py "RT_vel_Nrotate*_+z_*.dat" -ipix=256                   | (assuming channels = 101)

import argparse
import glob
import os
import numpy as np
from astropy.io import fits
from typing import Optional, List

def clamp_index(i: int, n: int) -> int:
    """Clamp i into [0, n-1]."""
    if i < 0:
        return 0
    if i >= n:
        return n - 1
    return i


def dat_to_fits(infile: str, ipix: int, channels_arg: Optional[int]) -> str:
    # Load columns
    data = np.loadtxt(infile)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < 4:
        raise ValueError(f"{infile}: expected >=4 columns (v, x, y, tr), got {data.shape[1]}")

    # csize = first row of 2nd column + last row of 2nd column
    x_first = float(data[0, 1])
    x_last = float(data[-1, 1])
    csize = x_first + x_last
    if csize == 0.0:
        raise ValueError(f"{infile}: computed csize=0 from x_first+x_last ({x_first}+{x_last})")

    # Velocity grid
    velocities = np.unique(data[:, 0])
    nvel = len(velocities)

    # Velocity channels
    channels = channels_arg if channels_arg is not None else nvel
    if nvel > channels:
        raise ValueError(
            f"{infile}: file contains {nvel} unique velocities but -channels={channels}. "
            f"Increase -channels or omit it."
        )

    # Allocate cube
    image_data = np.zeros((channels, ipix, ipix), dtype=np.float32)

    # Faster lookup: velocity value -> index
    v_to_idx = {v: i for i, v in enumerate(velocities)}

    # Fill
    scale = float(ipix) / float(csize)
    for v, x, y, tr in data[:, :4]:
        v_index = v_to_idx[float(v)]
        x_index = clamp_index(int(x * scale), ipix)
        y_index = clamp_index(int(y * scale), ipix)
        image_data[v_index, x_index, y_index] = tr

    # Output name: same as input but suffix .fits
    base, _ext = os.path.splitext(infile)
    outname = base + ".fits"

    # Header (minimal)
    header = fits.Header()
    header["COMMENT"] = "FITS cube created from (v, x, y, Tr) ASCII data"
    header["CUBESIZE"] = (csize, "Computed as x_first + x_last from column 2")
    header["IPIX"] = (ipix, "Spatial resolution (pixels)")
    header["NCHAN"] = (channels, "Velocity channels in output cube")

    fits.PrimaryHDU(data=image_data, header=header).writeto(outname, overwrite=True)
    return outname

def expand_patterns(patterns: List[str]) -> List[str]:
    files: list[str] = []
    for p in patterns:
        matches = sorted(glob.glob(p))
        if matches:
            files.extend(matches)
        else:
            files.append(p)
    seen = set()
    out = []
    for f in files:
        if f not in seen:
            seen.add(f)
            out.append(f)

    #omits the _cds.dat outputs
    filtered = [f for f in out if not f.endswith("_cds.dat")]

    return filtered


def main():
    ap = argparse.ArgumentParser(
        description="Convert RT_vel_* ASCII files (v, x, y, Tr) to FITS cubes."
    )
    ap.add_argument(
        "inputs",
        nargs="+",
        help="Input file(s) or glob pattern(s), e.g. RT_vel_Nrotate*.dat",
    )
    # -ipix=256, -channels=101
    ap.add_argument("-ipix", "--ipix", type=int, required=True, help="Spatial resolution (e.g. 256)")
    ap.add_argument(
        "-channels",
        "--channels",
        type=int,
        default=None,
        help="Number of velocity channels (e.g. 101). If omitted, inferred from file.",
    )

    args = ap.parse_args()
    files = expand_patterns(args.inputs)

    for infile in files:
        outname = dat_to_fits(infile, ipix=args.ipix, channels_arg=args.channels)
        print(f"{infile}  ->  {outname}")


if __name__ == "__main__":
    main()

