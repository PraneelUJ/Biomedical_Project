"""Export PICMUS simulated contrast/speckle UFF data in MATLAB-friendly HDF5."""
from pathlib import Path

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "PICMUS_simulation_contrast_speckle.uff"
TARGET = ROOT / "PICMUS_simulation_contrast_speckle_matlab.h5"

with h5py.File(SOURCE, "r") as src, h5py.File(TARGET, "w") as dst:
    channel = src["channel_data"]
    # UFF data is [transmit, receive, time]. MATLAB's h5read reverses HDF5
    # dimensions, so store [receive, time, transmit] for MATLAB [tx,time,rx].
    rf = channel["data"][()]
    dst.create_dataset("rf", data=np.transpose(rf, (1, 2, 0)),
                       compression="gzip", compression_opts=3)
    dst.create_dataset("sampling_frequency", data=channel["sampling_frequency"][()])
    dst.create_dataset("initial_time", data=channel["initial_time"][()])
    dst.create_dataset("sound_speed", data=channel["sound_speed"][()])
    geometry = channel["probe/geometry"][()]
    dst.create_dataset("element_x", data=geometry[0, :])
    dst.create_dataset("element_z", data=geometry[2, :])
    angles = np.array([
        src[f"channel_data/sequence/sequence_{i:04d}/source/azimuth"][()][0, 0]
        for i in range(1, rf.shape[0] + 1)
    ], dtype=np.float32)
    dst.create_dataset("angles", data=angles)
    dst.create_dataset("x_axis", data=src["scan/x_axis"][()].reshape(-1))
    dst.create_dataset("z_axis", data=src["scan/z_axis"][()].reshape(-1))
    nx = src["scan/x_axis"].shape[-1]
    nz = src["scan/z_axis"].shape[-1]
    ref = (src["beamformed_data/data/real"][()].reshape(-1) +
           1j * src["beamformed_data/data/imag"][()].reshape(-1))
    ref = ref.reshape((nz, nx), order="F")
    # Store transposed so MATLAB h5read returns [z,x].
    dst.create_dataset("reference_real", data=ref.real.T)
    dst.create_dataset("reference_imag", data=ref.imag.T)
    dst.attrs["source_file"] = SOURCE.name
    dst.attrs["source_md5"] = "d5d695a73ba0ee84db43a3281ddc7048"
    dst.attrs["raw_layout_in_matlab"] = "transmit,time,receive"

print(f"Exported {rf.shape[0]} tx x {rf.shape[2]} samples x {rf.shape[1]} rx to {TARGET}")
