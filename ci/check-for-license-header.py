"""
Check that all .py files in melodies_monet have the license header.
"""
from pathlib import Path

here = Path(__file__).parent
mm_lib = (here / "../melodies_monet").resolve()

header = """\
# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#
"""
header_lines = header.splitlines()

nbad = 0
for p in mm_lib.glob("**/*.py"):
    prel = p.relative_to(mm_lib)
    ok = True
    with open(p, "r") as f:
        for line, header_line in zip(f, header_lines):
            line = line.rstrip()
            if line != header_line:
                ok = False
                break

    if ok:
        print(prel, "ok")
    else:
        print(prel, "doesn't have correct header")
        print(f"  {line[:20]} ... != {header_line[:20]} ...")
        nbad += 1

raise SystemExit(1 if nbad > 0 else 0)
