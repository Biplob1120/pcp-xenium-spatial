"""Verify downloaded files against GEO filelist.txt sizes and basic readability.
GEO does not publish md5 sums for these supplementary files, so we check size + gzip/HDF5 integrity."""
import argparse
import gzip
import hashlib
from pathlib import Path

import h5py
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--data", required=True)
ap.add_argument("--md5", action="store_true", help="also write md5sums.txt for your own records")
a = ap.parse_args()
x = Path(a.data) / "xenium"
fl = pd.read_csv(x / "filelist.txt", sep="\t")
fl = fl[fl["#Archive/File"] == "File"].set_index("Name")
bad = []
for f in sorted(x.glob("GSM*")):
    if f.name not in fl.index:
        continue
    exp, got = int(fl.loc[f.name, "Size"]), f.stat().st_size
    ok = exp == got
    if ok and f.suffix == ".h5":
        try:
            h5py.File(f, "r")["matrix"]["shape"][:]
        except Exception as e:  # noqa
            ok = False
    elif ok and f.suffix == ".gz" and got < 3e9:
        try:
            with gzip.open(f) as fh:
                while fh.read(64 << 20):
                    pass
        except Exception:
            ok = False
    print(("OK   " if ok else "FAIL ") + f.name, got, exp)
    if not ok:
        bad.append(f.name)
    if a.md5:
        h = hashlib.md5(open(f, "rb").read()).hexdigest()
        with open(x / "md5sums.txt", "a") as out:
            out.write(f"{h}  {f.name}\n")
for t in sorted((Path(a.data) / "scrna").glob("*.tar.gz")):
    print("scRNA tar", t.name, t.stat().st_size)
print("FAILED:", bad if bad else "none")
raise SystemExit(1 if bad else 0)
