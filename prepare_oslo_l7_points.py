"""Export one L7 CIRS point-target frame and its scan metadata for MATLAB."""
from pathlib import Path

import h5py
import hdf5plugin  # register Open-H-RF's Blosc filter


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "oslo/C_verasonics_phantom/L7_FI_Verasonics_CIRS_points.hdf5"
RAW_OUT = ROOT / "oslo/C_verasonics_phantom/L7_FI_Verasonics_CIRS_points_frame1.h5"
META_OUT = ROOT / "oslo/C_verasonics_phantom/L7_FI_Verasonics_CIRS_points_metadata.h5"
METADATA_PATHS = {
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

with h5py.File(SOURCE, "r") as src:
    # Source layout: frame, transmit, time, receive element, channel.
    raw = src["/tracks/track_0/data/raw_data"][0, ..., 0]
    with h5py.File(RAW_OUT, "w") as dst:
        dst.attrs["source_file"] = SOURCE.name
        dst.attrs["frame_index"] = 0
        # h5read reverses HDF5 dimension order in MATLAB. Store reversed axes
        # so MATLAB receives the intended transmit x time x receive layout.
        dst.create_dataset("raw", data=raw.transpose(2, 1, 0),
                           compression="gzip", compression_opts=4)
    with h5py.File(META_OUT, "w") as dst:
        dst.attrs["source_file"] = SOURCE.name
        for name, path in METADATA_PATHS.items():
            value = src[path][()]
            # Compensate for MATLAB h5read's axis reversal, including the
            # square transmit-by-element delay matrix.
            if name in {"t0_delays", "tx_apodizations"}:
                value = value.T
            dst.create_dataset(name, data=value)

print(f"Exported raw frame {raw.shape} to {RAW_OUT}")
print(f"Exported scan metadata to {META_OUT}")
