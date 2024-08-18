#!/usr/bin/env python3

import pandas as pd
import sys
from lib.util import info

if __name__ == "__main__":
    sample_names = sys.argv[1].split(",")
    inpath = sys.argv[2]
    outpath = sys.argv[3]

    info(f"Loading depths from {inpath} for {len(sample_names)} samples into NetCDF: {outpath}")
    unitig_depth = pd.read_table(
        inpath, names=["unitig", *sample_names], index_col=["unitig"]
    ).rename_axis(columns="sample")
    unitig_depth.stack().to_xarray().to_netcdf(outpath)
