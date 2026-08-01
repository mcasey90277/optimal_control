#!/usr/bin/env python3
"""Generate a signature catalog of the pumpkyn / pumpkynPie MATLAB packages.

WHY THIS IS GENERATED, NOT WRITTEN. Both repos are external and move (two pulls
on 2026-07-30 alone). A hand-maintained catalogue of ~90 routines would drift,
and a doc that claims a signature which has since changed is worse than no doc
because it invites trust. This walks the packages and reads signatures and
purpose lines straight from the sources, so it can be regenerated and diffed.

The judgement layer -- what is relevant to our campaigns, what we have been
duplicating, what conventions bite -- lives in the hand-written companion
pumpkyn_reference.md, which this never touches.

Usage:
    python3 gen_pumpkyn_catalog.py > pumpkyn_catalog.md
"""
import re
import subprocess
import sys
from pathlib import Path

ROOTS = [
    ("pumpkyn", Path.home() / "Desktop/proj7/external/pumpkyn/src"),
    ("pumpkynPie", Path.home() / "Desktop/proj7/external/pumpkynPie/src"),
]

SIG = re.compile(r"^\s*function\s+(.*)$")


def head_commit(path: Path) -> str:
    """Short SHA of the repo containing `path`, or '?' if unavailable."""
    try:
        return subprocess.run(
            ["git", "-C", str(path), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=10,
        ).stdout.strip() or "?"
    except Exception:
        return "?"


def purpose(lines):
    """First substantive comment line after the signature.

    pumpkyn headers vary: some open with '%% Purpose:' then blank comment lines
    before the prose, others put the prose immediately after the signature.
    Skip section markers and empty comments and take the first real sentence.
    """
    for ln in lines:
        s = ln.strip()
        if not s.startswith("%"):
            if s and not s.startswith("function"):
                break          # into code; no header prose
            continue
        body = s.lstrip("%").strip()
        if not body:
            continue
        if body.lower().rstrip(":") in {"purpose", "inputs", "outputs", "notes"}:
            continue
        if set(body) <= {"-", "=", "_", "*"}:
            continue
        return body
    return ""


def main():
    print("# pumpkyn / pumpkynPie signature catalog")
    print()
    print("**GENERATED — do not edit.** Regenerate with")
    print("`python3 gen_pumpkyn_catalog.py > pumpkyn_catalog.md`.")
    print("Judgement and campaign relevance live in `pumpkyn_reference.md`.")
    print()
    for name, src in ROOTS:
        print(f"- `{name}` @ **{head_commit(src.parent)}**  (`{src.parent}`)")
    print()

    for name, src in ROOTS:
        if not src.exists():
            print(f"## {name}\n\n_not found at {src}_\n")
            continue
        print(f"## {name}")
        pkgs = sorted(p for p in src.rglob("+*") if p.is_dir())
        for pkg in pkgs:
            mfiles = sorted(pkg.glob("*.m"))
            if not mfiles:
                continue
            rel = "/".join(pkg.relative_to(src).parts)
            print(f"\n### `{rel}` — {len(mfiles)} routines\n")
            print("| routine | signature | purpose |")
            print("|---|---|---|")
            for m in mfiles:
                try:
                    lines = m.read_text(errors="replace").splitlines()
                except Exception:
                    continue
                sig, rest = "", []
                for i, ln in enumerate(lines[:40]):
                    hit = SIG.match(ln)
                    if hit:
                        sig = hit.group(1).strip().rstrip(";")
                        rest = lines[i + 1:i + 40]
                        break
                if not sig:
                    continue
                p = purpose(rest).replace("|", "\\|")
                if len(p) > 150:
                    p = p[:147] + "..."
                print(f"| `{m.stem}` | `{sig}` | {p} |")
        print()


if __name__ == "__main__":
    sys.exit(main())
