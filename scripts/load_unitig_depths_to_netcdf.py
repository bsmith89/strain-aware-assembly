#!/usr/bin/env python3

import xarray as xr
import sys
from lib.util import info
from tqdm import tqdm
import numpy as np

if __name__ == "__main__":
    sample_names = sys.argv[1].split(",")
    inpath = sys.argv[2]
    outpath = sys.argv[3]

    info(f"Loading depths from {inpath} for {len(sample_names)} samples into NetCDF: {outpath}")

    unitig_depth = []
    unitig_list = []
    with open(inpath) as f:
        for line in tqdm(f):
            unitig_id, *values = line.split()
            unitig_list.append(int(unitig_id))
            unitig_depth.append(np.array(values).astype(float))
    unitig_depth = xr.DataArray(np.stack(unitig_depth), coords=dict(unitig=unitig_list, sample=sample_names))

    unitig_depth.to_netcdf(outpath)
