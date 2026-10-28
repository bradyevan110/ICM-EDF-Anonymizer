# ICM EDF Anonymizer

Batch de-identification of EDF/EDF+ EEG recordings. Each file is renamed to its subject ID and
stripped of identifying information. The signal data is copied byte-for-byte, and the originals
are never modified.

## Requirements

Python 3.8+ with the standard library only. No packages to install.

## Folder layout

```
<input>/                     e.g. D:/ICM_EDFs/EBDM
└── sub-EMOP0462/
    └── <any name>.edf
```

Produces:

```
EDF_anon/                    de-identified files, ready to share
└── sub-EMOP0462/sub-EMOP0462.edf
EDF_anon_review/             annotation review notes (not for sharing)
└── sub-EMOP0462/annotations_review.txt
```

The subject folder name becomes the patient name in the output file.

## Usage

```bash
python deidentify_edf.py -i D:/ICM_EDFs/EBDM                         # output: EDF_anon next to the input folder
python deidentify_edf.py -i D:/ICM_EDFs/EBDM -o D:/somewhere/EDF_anon
python deidentify_edf.py -i D:/ICM_EDFs/EBDM --only sub-EMOP0462     # single subject
```

If `-i` is omitted, the input is the folder containing the script.

| Option | Effect |
|---|---|
| `--only SUB [SUB ...]` | Process only these subject folders |
| `--drop-sex` | Also remove sex from the patient field |
| `--clear-time` | Also set the start time to 00:00:00 |
| `--overwrite` | Redo subjects that already have output |

## What gets changed

| Field | After |
|---|---|
| Patient (MRN, sex, birth date, name) | `X <sex> X <subject-id>` |
| Recording (date, admin code, technician, equipment) | `Startdate X X X X` |
| Start date | `01.01.85`, the EDF+ convention for a removed date |
| Start time | kept (or cleared with `--clear-time`) |
| Annotation text | Montage names, patient name, MRN, birth date, technician name and anything that looks like a date are replaced with X's |
| Filename | `<subject-id>.edf` |

Annotation text is replaced with X's of the same length, so record layout and annotation timing
are unchanged.

## Converting to BIDS-iEEG

`edf_to_bids.py` turns the de-identified `EDF_anon` folder into a BIDS-iEEG dataset. It refuses
any file whose header has not been de-identified.

```bash
python edf_to_bids.py -i D:/ICM_EDFs/EDF_anon                        # output: BIDS next to the input folder
python edf_to_bids.py -i D:/ICM_EDFs/EDF_anon --task EBDM --session ieeg1 --reference "referential to ..."
```

Each recording becomes `sub-<id>/ses-<session>/ieeg/sub-<id>_ses-<session>_task-<task>_ieeg.edf`.
Next to it the script writes `_ieeg.json`, `_channels.tsv` (types inferred from the labels:
depth contacts as SEEG, 10-20 names as EEG, EKG, TRIG, DC inputs as ADC) and `_events.tsv`
(from the EDF+ annotations). It also writes `dataset_description.json`, `participants.tsv`
and `README`. Electrode coordinates are not included because the EDFs contain none.

| Option | Default |
|---|---|
| `--task` | `EBDM` |
| `--session` | `ieeg1` |
| `--line-freq` | `60` |
| `--reference` | `n/a` (fill in the recording reference) |
| `--institution` | `n/a` |
| `--only`, `--overwrite` | as above |

## Before sharing

Free-text annotations can contain anything a clinician typed. Read each
`annotations_review.txt` before sharing; it lists every unique annotation text left in the file.
