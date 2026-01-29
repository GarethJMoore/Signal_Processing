# TAS data pipeline

`tas_pipeline.py` turns the notebook logic into a repeatable pipeline for experimenters who need to process new scans. It loads the raw scans, averages them, performs negative-time background subtraction, applies chirp correction, and writes out reusable artifacts for downstream analysis.

## Requirements

- Python 3.10+
- `numpy`
- `scipy`

## Usage

Run the pipeline against the provided scans:

```bash
python tas_pipeline.py --data-dir ../Data --output-dir ../outputs
```

Optional flags:

```bash
python tas_pipeline.py \
  --data-dir ../Data \
  --output-dir ../outputs \
  --background-time-threshold 0 \
  --save-csv
```

## Outputs

The pipeline writes `.npy` arrays and a `metadata.json` summary to the output directory. Use `--save-csv` to also emit CSV versions of the time, wavelength, and chirp-corrected data.
