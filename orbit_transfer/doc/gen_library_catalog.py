#!/usr/bin/env python3
"""Generate library_catalog.md — the function reference for OUR orbit-transfer
library layers (oclib/+oc, costate_common, verify_common, cr3bp_common).

Mirrors the gen_pumpkyn_catalog.py pattern: parse each .m file's signature
line and its pumpkyn-style header (%% Purpose block, or the H1 comment line
for flat-style files), emit one markdown catalog. Regenerate after any
library change; the file carries its own regeneration command so it cannot
silently drift.

Usage:  python3 gen_library_catalog.py   (from anywhere; paths are absolute
        relative to this script's location)
"""
import re
from pathlib import Path
from datetime import date

HERE = Path(__file__).resolve().parent            # orbit_transfer/doc
OT = HERE.parent                                   # orbit_transfer
ROOT = OT.parent                                   # optimal_control

LAYERS = [
    ("oclib/+oc", ROOT / "oclib" / "+oc",
     "Cross-top-level-folder core (consumers in orbit_transfer AND "
     "booster_landing). Admission: second top-level consumer + equivalence "
     "gate. See ../../oclib/README.md for consumers/gates per function."),
    ("costate_common", OT / "costate_common",
     "The costate-pipeline library: family construction, multiple shooting, "
     "seeds, conjugate test, catalog build/validate/sweep. See its "
     "README.md for the judgment layer."),
    ("verify_common", OT / "verify_common",
     "First-order optimality gate layer + the shared continuous-residual "
     "(G1) gate. See its README.md and OPTIMALITY_CERTIFICATION.md."),
    ("cr3bp_common", OT / "cr3bp_common",
     "Shared CR3BP GTO problem definition (params, endpoints, setup)."),
]


def parse_mfile(path: Path):
    """Return (signature, purpose, inputs, outputs) from a MATLAB file."""
    text = path.read_text(errors="replace")
    lines = text.splitlines()
    sig = ""
    for ln in lines:
        s = ln.strip()
        if s.startswith("function"):
            sig = re.sub(r"^function\s+", "", s).rstrip(";")
            break
    # Purpose: pumpkyn '%% Purpose:' block, else first contiguous comment
    purpose = []
    m = re.search(r"%%\s*Purpose:\s*\n(.*?)(?=\n%%|\n[^%])", text, re.S)
    if m:
        block = m.group(1)
    else:
        block = ""
        started = False
        for ln in lines[1:]:
            if ln.strip().startswith("%"):
                block += ln + "\n"
                started = True
            elif started:
                break
    for ln in block.splitlines():
        t = ln.strip().lstrip("%").strip()
        if t:
            purpose.append(t)
        elif purpose:
            break  # first paragraph only
    # Inputs/Outputs names (pumpkyn column style or 'name - desc' style)
    def io_names(section):
        mm = re.search(rf"%%?\s*(?:{section}):?\s*\n(.*?)(?=\n%%|\n\s*% *(?:OUTPUTS|REFERENCES|Revision)|\nfunction|\n[^%])",
                       text, re.S | re.I)
        if not mm:
            return []
        names = []
        for ln in mm.group(1).splitlines():
            t = ln.strip().lstrip("%").rstrip()
            m2 = re.match(r"^\s*\.?([A-Za-z]\w*)\s{2,}", t) or \
                 re.match(r"^\s*\.?([A-Za-z]\w*)\s+-\s", t)
            if m2:
                nm = m2.group(1)
                if nm not in names:
                    names.append(nm)
        return names
    return sig, " ".join(purpose), io_names("Inputs"), io_names("OUTPUTS|Outputs")


def main():
    out = [
        "# Orbit-transfer library catalog — generated function reference",
        "",
        f"Generated {date.today().isoformat()} by `gen_library_catalog.py`.",
        "**Do not edit by hand** — regenerate with:",
        "```sh",
        "python3 orbit_transfer/doc/gen_library_catalog.py",
        "```",
        "Per-function headers in the .m files remain the authoritative",
        "documentation (derivations, assumptions, revision history); this",
        "catalog is the browsable index for humans and Claude sessions.",
        "The judgment layer (what to use when, consumers, gates, traps)",
        "lives in each layer's README.md and in `../OCP_UNIFYING_MATH.md`.",
        "",
    ]
    for name, folder, blurb in LAYERS:
        out += [f"## {name}", "", blurb, ""]
        files = sorted(p for p in folder.glob("*.m"))
        if not files:
            out += ["*(no functions found)*", ""]
            continue
        for f in files:
            sig, purpose, ins, outs = parse_mfile(f)
            out.append(f"### `{f.name}`")
            if sig:
                out.append(f"`{sig}`  ")
            if purpose:
                out.append(purpose)
            io = []
            if ins:
                io.append("in: " + ", ".join(f"`{n}`" for n in ins))
            if outs:
                io.append("out: " + ", ".join(f"`{n}`" for n in outs))
            if io:
                out.append("*" + " · ".join(io) + "*")
            out.append("")
        # tests subfolder, names only
        tdir = folder / "tests"
        if tdir.is_dir():
            tnames = ", ".join(f"`{p.name}`" for p in sorted(tdir.glob("*.m")))
            if tnames:
                out += [f"**tests/**: {tnames}", ""]
    dst = HERE / "library_catalog.md"
    dst.write_text("\n".join(out) + "\n")
    n_funcs = sum(1 for L in out if L.startswith("### "))
    print(f"wrote {dst} ({n_funcs} functions, {len(out)} lines)")


if __name__ == "__main__":
    main()
