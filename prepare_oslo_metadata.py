"""Export compressed Oslo acquisition metadata to plain HDF5 for MATLAB h5read."""
from pathlib import Path
import h5py
import hdf5plugin  # registers Open-H-RF's Blosc filter (32001)

root = Path(__file__).resolve().parent
source = root / "oslo/C_verasonics_phantom/FI_P4_cysts_center.hdf5"
target = root / "oslo/C_verasonics_phantom/FI_P4_cysts_center_metadata.h5"
paths = {
    "sampling_frequency": "/tracks/track_0/scan/sampling_frequency",
    "center_frequency": "/tracks/track_0/scan/center_frequency",
    "sound_speed": "/tracks/track_0/scan/sound_speed",
    "initial_times": "/tracks/track_0/scan/initial_times",
    "t0_delays": "/tracks/track_0/scan/t0_delays",
    "focus_distances": "/tracks/track_0/scan/focus_distances",
    "polar_angles": "/tracks/track_0/scan/polar_angles",
    "transmit_origins": "/tracks/track_0/scan/transmit_origins",
    "tx_apodizations": "/tracks/track_0/scan/tx_apodizations",
    "probe_geometry": "/probe/probe_geometry",
}
with h5py.File(source, "r") as src, h5py.File(target, "w") as dst:
    dst.attrs["source_file"] = source.name
    for name, path in paths.items():
        dst.create_dataset(name, data=src[path][()])  # no compression/plugin
print(f"Wrote plain MATLAB-readable metadata: {target}")

# Also prepare one point-scatterer frame. Gzip is supported by MATLAB's built-in
# HDF5 reader, unlike the source file's Blosc filter.
points_source = root / "oslo/C_verasonics_phantom/FI_P4_point_scatterers.hdf5"
points_target = root / "oslo/C_verasonics_phantom/FI_P4_point_scatterers_frame1.h5"
points_meta = root / "oslo/C_verasonics_phantom/FI_P4_point_scatterers_metadata.h5"
with h5py.File(points_source, "r") as src:
    with h5py.File(points_target, "w") as dst:
        frame = src["/tracks/track_0/data/raw_data"][0, ..., 0]
        dst.create_dataset("raw", data=frame, compression="gzip", compression_opts=4)
    with h5py.File(points_meta, "w") as dst:
        dst.attrs["source_file"] = points_source.name
        for name, path in paths.items():
            dst.create_dataset(name, data=src[path][()])
print(f"Wrote point-scatterer frame and plain metadata: {points_target}")
