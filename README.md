# Biomedical ultrasound RF denoising

MATLAB experiments for ultrasound channel-RF simulation, delay-and-sum (DAS)
beamforming, and noise reduction. The code covers a Field II cyst phantom,
Open-H-RF Oslo acquisitions, and the PICMUS simulated contrast/speckle data.
Dataset-specific scripts use each acquisition's geometry and timing metadata;
their SNR results should be interpreted within that dataset and reference.

## Source code

| File | Purpose |
|---|---|
| `fieldii_cyst_sim.m` | Generate a Field II cyst/point-target RF simulation and save `fieldII_cyst_RF.mat`. |
| `DAS_adaptive_CF.m` | DAS, coherence-factor, and filter-bank comparisons on the Field II output. |
| `adaptive_frequency_bank.m` | Reusable overlapping-band adaptive RF filtering function. |
| `oslo_adaptive_das.m` | DAS and adaptive denoising experiment for Oslo P4 data. |
| `oslo_l7_points_experiment.m` | DAS and post-beamforming Wiener experiment on Oslo L7 CIRS point targets. |
| `PICMUS_contrast_denoising_experiment.m` | PICMUS DAS comparison with the original post-DAS Wiener method. |
| `PICMUS_filterbank_sweep.m` | Separate parameter sweep across RF banks, Wiener settings, and cascaded methods. |
| `prepare_oslo_metadata.py` | Export MATLAB-readable P4 metadata and a point-scatterer frame from Open-H-RF Blosc HDF5. |
| `prepare_oslo_l7_points.py` | Export MATLAB-readable L7 CIRS point-target RF and metadata. |
| `prepare_picmus_contrast.py` | Convert the PICMUS UFF acquisition to MATLAB-friendly HDF5. |

Generated `.mat`, `.h5`, `.uff`, and sweep CSV files are local outputs and are
ignored by Git. The selected Oslo `.hdf5` acquisitions already in the
repository use Git LFS. PICMUS data is not included; download the simulated
contrast/speckle UFF file from the [PICMUS dataset record on Zenodo](https://zenodo.org/records/20261898)
and place it beside the scripts as `PICMUS_simulation_contrast_speckle.uff`.

## Setup

Install the Python conversion dependencies:

```bash
python3 -m pip install -r requirements.txt
```

Field II is an external MATLAB toolbox and is not vendored here. Install it
separately, then add its `m_files` directory to the MATLAB path before running
the Field II scripts. MATLAB's Signal Processing Toolbox is also needed for
analytic-signal and resampling functions (`hilbert`, `resample`).

## Run experiments

Run commands from the project root. For Oslo P4 metadata and a point-scatterer
frame, first run:

```bash
python3 prepare_oslo_metadata.py
```

Then run `oslo_adaptive_das.m` in MATLAB. For the L7 CIRS point-target
acquisition, run `python3 prepare_oslo_l7_points.py`, then run
`oslo_l7_points_experiment.m`.

For PICMUS, after placing the UFF data in the project root, run:

```bash
python3 prepare_picmus_contrast.py
```

Then run `PICMUS_contrast_denoising_experiment.m` for the Wiener-only result or
`PICMUS_filterbank_sweep.m` for the longer comparison. The sweep evaluates
multiple passbands, 2/3/4/6/8-band RF banks, local Wiener settings, RF-bank plus
post-DAS Wiener cascades, and four repeatable noise realizations. It saves a
ranked CSV and MATLAB results file after completion.

For Field II, run `fieldii_cyst_sim.m` first, then `DAS_adaptive_CF.m`.

## Interpretation

The experiments add controlled Gaussian noise to channel RF and evaluate
reconstructed images against a clean reconstruction using matched processing.
The PICMUS sweep ranks methods by error SNR and also records envelope
correlation and mean envelope ratio, since aggressive suppression can improve
error SNR while reducing image amplitude. These synthetic-noise experiments
do not establish performance on scanner noise.
