"""
Batch de-identification of EDF/EDF+ files.

Expected layout:   <input_dir>/<subject-id>/<anything>.edf
Output layout:     <output_dir>/<subject-id>/<subject-id>.edf

For each file it:
  * Replaces the patient field with  "X <sex> X <subject-id>"  (MRN, DOB and name removed;
    sex kept unless --drop-sex).
  * Replaces the recording field with "Startdate X X X X" (removes recording date,
    admin code, technician and equipment).
  * Sets the header start date to 01.01.85 (the EDF+ convention for an anonymized date).
    Start time of day is kept unless --clear-time.
  * Scrubs free-text in the "EDF Annotations" channel: montage names, patient name tokens, MRN, DOB,
    technician name and anything that looks like a calendar date are overwritten with X's
    (same length, so record layout and annotation timing are unchanged).
  * Renames the file to the subject ID (original filenames contain patient names).
  * Writes <output_dir>_review/<subject-id>/annotations_review.txt listing every unique
    annotation text left in the file, so a human can confirm nothing identifying remains.
    (Kept outside <output_dir> so that folder holds only the EDFs to share.)

Signal data are copied byte-for-byte. Originals are never modified.

Usage:
    python deidentify_edf.py                       # input = folder containing this script,
                                                   # output = EDF_anon next to it
    python deidentify_edf.py -i D:/ICM_EDFs/EBDM -o D:/ICM_EDFs/EDF_anon
    python deidentify_edf.py --only sub-EMOP0462   # run a single subject
"""
import argparse
import collections
import re
import sys
from pathlib import Path

CHUNK_RECORDS = 2000  # data records processed per read/write

DATE_PATTERNS = [
    r"\b\d{1,4}[-/.]\d{1,2}[-/.]\d{2,4}\b",                    # 5-4-2026, 05/04/26, 2026-05-04
    r"\b\d{1,2}[- ]?(?:JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)[A-Z]*[- ,]*\d{2,4}\b",
    r"\b(?:JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)[A-Z]*\.? \d{1,2},? \d{2,4}\b",
]


def field(h, start, length):
    return h[start:start + length].decode("latin-1")


def pad(text, length):
    b = text.encode("ascii", "replace")
    if len(b) > length:
        raise ValueError(f"'{text}' is longer than {length} bytes")
    return b.ljust(length, b" ")


def phi_tokens(patient, recording, src_name):
    """Collect identifying strings from the original header and filename."""
    tokens = set()
    p = patient.split(" ")
    if len(p) >= 4:  # EDF+: code sex birthdate name
        code, birthdate, name = p[0], p[2], " ".join(p[3:])
        tokens.update(t for t in (code, birthdate) if t != "X")
    else:
        name = patient
    tokens.update(re.split(r"[\s,_^~]+", name))
    r = recording.split(" ")
    if len(r) >= 5 and r[0] == "Startdate":  # Startdate date admincode technician equipment
        tokens.update(re.split(r"[\s,_^~]+", r[3]))
    stem = src_name.rsplit("_", 1)[0]  # "Last~ First_<uuid>.EDF" -> "Last~ First"
    tokens.update(re.split(r"[\s,_^~]+", stem))
    return {t for t in tokens if len(t) >= 2 and t.upper() not in {"X", "M", "F"}}


def build_scrubber(tokens):
    # short tokens only as whole words so e.g. "Al" doesn't hit "Alpha"
    parts = [re.escape(t) if len(t) >= 4 else rf"\b{re.escape(t)}\b"
             for t in sorted(tokens, key=len, reverse=True)]
    # Natus montage names embed the patient name, sometimes spelled differently from the header,
    # so blank the whole name (first line after "Montage:")
    montage = r"(?<=Montage:)[^\r\n]+"
    regex = re.compile("|".join([montage] + parts + DATE_PATTERNS), re.IGNORECASE)
    return lambda s: regex.subn(lambda m: "X" * len(m.group()), s)


def scrub_annotations(raw, scrub, seen):
    """Scrub TAL text fields in one record's annotation bytes, keeping length identical."""
    out, n = [], 0
    for tal in raw.split(b"\x00"):
        if b"\x14" in tal:
            onset, *texts = tal.split(b"\x14")
            new = []
            for t in texts:
                s, k = scrub(t.decode("latin-1"))
                n += k
                if s:
                    seen[s] += 1
                new.append(s.encode("latin-1"))
            tal = b"\x14".join([onset] + new)
        out.append(tal)
    result = b"\x00".join(out)
    assert len(result) == len(raw)
    return result, n


def process(src, dst, review, subject, drop_sex, clear_time):
    with open(src, "rb") as fin:
        hdr = bytearray(fin.read(256))
        ns = int(field(hdr, 252, 4))
        hdr += fin.read(ns * 256)
        if len(hdr) != int(field(hdr, 184, 8)):
            raise ValueError("unexpected header size")
        nrec = int(field(hdr, 236, 8))

        patient, recording = field(hdr, 8, 80).strip(), field(hdr, 88, 80).strip()
        tokens = phi_tokens(patient, recording, src.name)
        scrub = build_scrubber(tokens)

        p = patient.split(" ")
        sex = p[1] if len(p) >= 4 and p[1] in ("M", "F") and not drop_sex else "X"
        hdr[8:88] = pad(f"X {sex} X {subject}", 80)
        hdr[88:168] = pad("Startdate X X X X", 80)
        hdr[168:176] = b"01.01.85"
        if clear_time:
            hdr[176:184] = b"00.00.00"

        sig = lambda i, off, size: field(hdr, 256 + ns * off + i * size, size)
        labels = [sig(i, 0, 16).strip() for i in range(ns)]
        nsamp = [int(sig(i, 216, 8)) for i in range(ns)]
        recsize = sum(nsamp) * 2
        annots = [(sum(nsamp[:i]) * 2, nsamp[i] * 2) for i, l in enumerate(labels) if l == "EDF Annotations"]

        # check the new patient/recording fields and all signal headers (not the date/time fields)
        leftover = scrub(field(hdr, 8, 160) + hdr[256:].decode("latin-1"))[1]
        if leftover:
            raise RuntimeError("identifying text found in signal headers; inspect this file manually")

        seen = collections.Counter()
        redactions = 0
        dst.parent.mkdir(parents=True, exist_ok=True)
        tmp = dst.with_suffix(".partial")
        with open(tmp, "wb") as fout:
            fout.write(hdr)
            done = 0
            while done < nrec:
                n = min(CHUNK_RECORDS, nrec - done)
                buf = bytearray(fin.read(n * recsize))
                if len(buf) != n * recsize:
                    raise ValueError("file is shorter than its header claims")
                for r in range(n):
                    for off, size in annots:
                        a = r * recsize + off
                        buf[a:a + size], k = scrub_annotations(bytes(buf[a:a + size]), scrub, seen)
                        redactions += k
                fout.write(buf)
                done += n
                print(f"\r  {subject}: {100 * done / nrec:5.1f}%", end="", flush=True)
        tmp.replace(dst)
    print()

    review.parent.mkdir(parents=True, exist_ok=True)
    with open(review, "w", encoding="utf-8") as f:
        f.write(f"{subject}: {redactions} redaction(s) in annotations\n")
        f.write("Unique annotation texts after scrubbing (count, text):\n\n")
        for text, c in seen.most_common():
            f.write(f"{c:8d}  {text!r}\n")
    return redactions, len(seen)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    here = Path(__file__).resolve().parent
    ap.add_argument("-i", "--input", type=Path, default=here)
    ap.add_argument("-o", "--output", type=Path, default=None, help="default: EDF_anon next to the input folder")
    ap.add_argument("--only", nargs="*", help="subject folder name(s) to process")
    ap.add_argument("--drop-sex", action="store_true", help="also remove sex from the patient field")
    ap.add_argument("--clear-time", action="store_true", help="also set the start time to 00.00.00")
    ap.add_argument("--overwrite", action="store_true", help="redo subjects that already have output")
    args = ap.parse_args()
    out_root = args.output or args.input.resolve().parent / "EDF_anon"
    review_root = out_root.with_name(out_root.name + "_review")

    jobs = []
    for sub_dir in sorted(d for d in args.input.iterdir() if d.is_dir()):
        if args.only and sub_dir.name not in args.only:
            continue
        files = sorted(f for f in sub_dir.iterdir() if f.suffix.lower() == ".edf")
        for k, f in enumerate(files, 1):
            suffix = f"_run-{k:02d}" if len(files) > 1 else ""
            jobs.append((f, out_root / sub_dir.name / f"{sub_dir.name}{suffix}.edf",
                         review_root / sub_dir.name / f"annotations_review{suffix}.txt", sub_dir.name))

    print(f"{len(jobs)} file(s) -> {out_root}")
    failed = []
    for src, dst, review, subject in jobs:
        if dst.exists() and not args.overwrite:
            print(f"  {subject}: output exists, skipping (use --overwrite to redo)")
            continue
        try:
            red, uniq = process(src, dst, review, subject, args.drop_sex, args.clear_time)
            print(f"  {subject}: done, {red} annotation redaction(s), {uniq} unique annotation text(s)")
        except Exception as e:
            print(f"\n  {subject}: FAILED - {e}")
            failed.append(subject)
    if failed:
        print("Failed:", ", ".join(failed))
        sys.exit(1)
    print(f"All done. Review the files in {review_root} before sharing.")


if __name__ == "__main__":
    main()
