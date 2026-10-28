"""
Convert de-identified EDF+ recordings into a BIDS-iEEG dataset.

Expected input (output of deidentify_edf.py):
    <input_dir>/sub-<label>/sub-<label>.edf             (or sub-<label>_run-01.edf, ...)

Output:
    <output_dir>/
        dataset_description.json, README, participants.tsv, participants.json
        sub-<label>/ses-<session>/ieeg/
            sub-<label>_ses-<session>_task-<task>[_run-01]_ieeg.edf
            ..._ieeg.json       recording metadata
            ..._channels.tsv    one row per signal (type, units, sampling rate)
            ..._events.tsv      EDF+ annotations (onset in seconds from recording start)
            ..._events.json

Channel types are inferred from the labels (see channel_type()). Depth contacts such as
"RAm1" become SEEG; check the printed warnings for any label that falls through to MISC.

Electrode coordinates (_electrodes.tsv / _coordsystem.json) are not written: the EDF files
contain no positions. Add them from imaging if available.

The input must already be de-identified: files whose patient field doesn't name the subject
folder, or whose recording field still has content, are refused.

Signal data are copied byte-for-byte. Input files are never modified.

Usage:
    python edf_to_bids.py -i D:/ICM_EDFs/EDF_anon                    # output: BIDS next to the input
    python edf_to_bids.py -i D:/ICM_EDFs/EDF_anon -o D:/ICM_EDFs/BIDS --task EBDM --session ieeg1
"""
import argparse
import json
import re
import sys
from pathlib import Path

BIDS_VERSION = "1.10.0"
CHUNK_RECORDS = 2000  # data records processed per read/write

SCALP_EEG = {"FP1", "FP2", "FPZ", "F3", "F4", "F7", "F8", "FZ", "C3", "C4", "CZ", "P3", "P4",
             "P7", "P8", "PZ", "O1", "O2", "OZ", "T3", "T4", "T5", "T6", "T7", "T8", "A1", "A2",
             "M1", "M2", "F9", "F10", "T9", "T10", "P9", "P10"}


def channel_type(label):
    """Return (BIDS type, description) for a Natus channel label."""
    up = label.upper()
    if up in SCALP_EEG:
        return "EEG", "scalp EEG"
    if up in ("EKG", "ECG") or re.fullmatch(r"(EKG|ECG)\d*", up):
        return "ECG", "electrocardiogram"
    if re.fullmatch(r"E[MO]G\d*", up):
        return up[:3], ""
    if up == "TRIG":
        return "TRIG", "trigger channel"
    if re.fullmatch(r"DC\d+", up):
        return "ADC", "DC analog input"
    if up == "PLETH":
        return "MISC", "plethysmograph"
    if up == "OSAT":
        return "MISC", "oxygen saturation"
    if up == "PR":
        return "MISC", "pulse rate"
    if re.fullmatch(r"C\d{3}", up):
        return "MISC", "unlabeled amplifier input"
    if re.fullmatch(r"[LR][A-Za-z]+\d+", label):
        return "SEEG", ""
    return "MISC", ""


def field(h, start, length):
    return h[start:start + length].decode("latin-1")


def read_header(path):
    with open(path, "rb") as f:
        hdr = f.read(256)
        ns = int(field(hdr, 252, 4))
        hdr += f.read(ns * 256)
    sig = lambda off, size: [field(hdr, 256 + ns * off + i * size, size).strip() for i in range(ns)]
    return {
        "raw": hdr,
        "patient": field(hdr, 8, 80).strip(),
        "recording": field(hdr, 88, 80).strip(),
        "header_bytes": int(field(hdr, 184, 8)),
        "nrec": int(field(hdr, 236, 8)),
        "rec_dur": float(field(hdr, 244, 8)),
        "labels": sig(0, 16),
        "units": sig(96, 8),
        "prefilter": sig(136, 80),
        "nsamp": [int(x) for x in sig(216, 8)],
    }


def check_deidentified(h, subject_dir):
    p = h["patient"].split(" ")
    if p[-1] != subject_dir or h["recording"] != "Startdate X X X X":
        raise RuntimeError("header does not look de-identified; run deidentify_edf.py first")


def parse_tals(raw):
    """Yield (onset, duration, [texts]) for each TAL in one record's annotation bytes."""
    for tal in raw.split(b"\x00"):
        if b"\x14" not in tal:
            continue
        timing, *texts = tal.decode("latin-1").split("\x14")
        onset, _, duration = timing.partition("\x15")
        yield float(onset), float(duration) if duration else 0.0, [t for t in texts if t]


def copy_and_extract_events(src, dst, h):
    """Copy the EDF to dst and return its annotations as (onset, duration, text) tuples."""
    recsize = sum(h["nsamp"]) * 2
    annots = [(sum(h["nsamp"][:i]) * 2, h["nsamp"][i] * 2)
              for i, l in enumerate(h["labels"]) if l == "EDF Annotations"]
    events, start = [], None
    tmp = dst.with_suffix(".partial")
    with open(src, "rb") as fin, open(tmp, "wb") as fout:
        fout.write(fin.read(h["header_bytes"]))
        done = 0
        while done < h["nrec"]:
            n = min(CHUNK_RECORDS, h["nrec"] - done)
            buf = fin.read(n * recsize)
            if len(buf) != n * recsize:
                raise ValueError("file is shorter than its header claims")
            for r in range(n):
                for off, size in annots:
                    a = r * recsize + off
                    for k, (onset, duration, texts) in enumerate(parse_tals(buf[a:a + size])):
                        if start is None and k == 0:
                            start = onset  # time-keeping TAL of the first record
                        events.extend((onset, duration, t) for t in texts)
            fout.write(buf)
            done += n
            print(f"\r  {dst.name}: {100 * done / h['nrec']:5.1f}%", end="", flush=True)
    tmp.replace(dst)
    print()
    start = start or 0.0
    return sorted((o - start, d, t) for o, d, t in events)


def write_tsv(path, columns, rows):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(columns) + "\n")
        for row in rows:
            f.write("\t".join(str(row[c]) for c in columns) + "\n")


def write_json(path, obj):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)
        f.write("\n")


def fmt(x):
    return f"{x:.6f}".rstrip("0").rstrip(".")


def write_sidecars(stem, h, events, args):
    data = [i for i, l in enumerate(h["labels"]) if l != "EDF Annotations"]
    fs = {i: h["nsamp"][i] / h["rec_dur"] for i in data}
    main_fs = max(set(fs.values()), key=list(fs.values()).count)

    rows, counts = [], {}
    for i in data:
        label = h["labels"][i]
        ctype, desc = channel_type(label)
        if ctype == "MISC" and not desc:
            print(f"    warning: channel '{label}' not recognised, typed as MISC")
        counts[ctype] = counts.get(ctype, 0) + 1
        hp = re.search(r"HP:\s*([\d.]+)", h["prefilter"][i])
        lp = re.search(r"LP:\s*([\d.]+)", h["prefilter"][i])
        rows.append({
            "name": label,
            "type": ctype,
            "units": h["units"][i] or "n/a",
            "low_cutoff": hp.group(1) if hp else "n/a",
            "high_cutoff": lp.group(1) if lp else "n/a",
            "sampling_frequency": fmt(fs[i]),
            "status": "n/a",
            "description": desc or "n/a",
        })
    write_tsv(f"{stem}_channels.tsv", list(rows[0]), rows)

    write_json(f"{stem}_ieeg.json", {
        "TaskName": args.task,
        "InstitutionName": args.institution,
        "Manufacturer": "Natus",
        "PowerLineFrequency": args.line_freq,
        "SamplingFrequency": main_fs,
        "SoftwareFilters": "n/a",
        "iEEGReference": args.reference,
        "RecordingType": "continuous",
        "RecordingDuration": round(h["nrec"] * h["rec_dur"], 6),
        "SEEGChannelCount": counts.get("SEEG", 0),
        "ECOGChannelCount": counts.get("ECOG", 0),
        "EEGChannelCount": counts.get("EEG", 0),
        "ECGChannelCount": counts.get("ECG", 0),
        "EMGChannelCount": counts.get("EMG", 0),
        "EOGChannelCount": counts.get("EOG", 0),
        "TriggerChannelCount": counts.get("TRIG", 0),
        "MiscChannelCount": counts.get("MISC", 0) + counts.get("ADC", 0),
    })

    ev_rows = []
    for onset, duration, text in events:
        if text.startswith("Montage:"):  # montage names were blanked during de-identification
            continue
        text = re.sub(r"[\t\r\n]+", " ", text).strip()
        trig = re.fullmatch(r"TRIG\[(\d+)\]:?", text)
        ev_rows.append({"onset": fmt(onset), "duration": fmt(duration),
                        "trial_type": text or "n/a", "value": int(trig.group(1)) if trig else "n/a"})
    write_tsv(f"{stem}_events.tsv", ["onset", "duration", "trial_type", "value"], ev_rows)
    write_json(f"{stem}_events.json", {
        "trial_type": {"Description": "Annotation text from the EDF+ annotation channel"},
        "value": {"Description": "Trigger code for TRIG[n] annotations"},
    })
    return len(ev_rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", type=Path, required=True, help="folder of de-identified subject folders")
    ap.add_argument("-o", "--output", type=Path, default=None, help="default: BIDS next to the input folder")
    ap.add_argument("--task", default="EBDM", help="BIDS task label (default: EBDM)")
    ap.add_argument("--session", default="ieeg1", help="BIDS session label (default: ieeg1)")
    ap.add_argument("--name", default="EBDM", help="dataset name for dataset_description.json")
    ap.add_argument("--line-freq", type=int, default=60, help="power line frequency in Hz (default: 60)")
    ap.add_argument("--reference", default="n/a", help="iEEGReference text, e.g. 'referential to RAm1'")
    ap.add_argument("--institution", default="n/a", help="InstitutionName for _ieeg.json")
    ap.add_argument("--only", nargs="*", help="subject folder name(s) to process")
    ap.add_argument("--overwrite", action="store_true", help="redo recordings already in the output")
    args = ap.parse_args()
    for label in (args.task, args.session):
        if not re.fullmatch(r"[A-Za-z0-9]+", label):
            sys.exit(f"BIDS labels must be letters and digits only: '{label}'")
    out_root = args.output or args.input.resolve().parent / "BIDS"
    out_root.mkdir(parents=True, exist_ok=True)

    participants, failed = {}, []
    for sub_dir in sorted(d for d in args.input.iterdir() if d.is_dir() and d.name.startswith("sub-")):
        label = re.sub(r"[^A-Za-z0-9]", "", sub_dir.name[4:])
        for src in sorted(f for f in sub_dir.iterdir() if f.suffix.lower() == ".edf"):
            try:
                h = read_header(src)
                check_deidentified(h, sub_dir.name)
                p = h["patient"].split(" ")
                participants[label] = p[1] if len(p) >= 4 and p[1] in ("M", "F") else "n/a"
                if args.only and sub_dir.name not in args.only:
                    continue
                run = re.search(r"_run-(\d+)", src.stem)
                ieeg_dir = out_root / f"sub-{label}" / f"ses-{args.session}" / "ieeg"
                ieeg_dir.mkdir(parents=True, exist_ok=True)
                stem = ieeg_dir / (f"sub-{label}_ses-{args.session}_task-{args.task}"
                                   + (f"_run-{run.group(1)}" if run else ""))
                dst = Path(f"{stem}_ieeg.edf")
                if dst.exists() and not args.overwrite:
                    print(f"  {dst.name}: exists, skipping (use --overwrite to redo)")
                    continue
                events = copy_and_extract_events(src, dst, h)
                n_ev = write_sidecars(stem, h, events, args)
                print(f"  sub-{label}: done, {n_ev} event(s)")
            except Exception as e:
                print(f"\n  {src.name}: FAILED - {e}")
                failed.append(src.name)

    write_json(out_root / "dataset_description.json", {
        "Name": args.name,
        "BIDSVersion": BIDS_VERSION,
        "DatasetType": "raw",
        "GeneratedBy": [{"Name": "ICM EDF Anonymizer (edf_to_bids.py)"}],
    })
    write_tsv(out_root / "participants.tsv", ["participant_id", "sex", "age"],
              [{"participant_id": f"sub-{k}", "sex": v, "age": "n/a"} for k, v in sorted(participants.items())])
    write_json(out_root / "participants.json", {
        "sex": {"Description": "sex of the participant", "Levels": {"M": "male", "F": "female"}},
        "age": {"Description": "age of the participant", "Units": "years"},
    })
    readme = out_root / "README"
    if not readme.exists():
        readme.write_text(f"{args.name}\n\nDe-identified intracranial EEG (sEEG) recordings in BIDS-iEEG format.\n"
                          "Recording dates were removed during de-identification (EDF start date set to 01.01.85).\n",
                          encoding="utf-8")
    if failed:
        print("Failed:", ", ".join(failed))
        sys.exit(1)
    print(f"All done -> {out_root}")


if __name__ == "__main__":
    main()
