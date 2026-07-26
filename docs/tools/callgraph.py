#!/usr/bin/env python3
"""callgraph.py -- static call-graph survey for a MATLAB campaign folder.

Usage:
    python3 docs/tools/callgraph.py <root-folder> [entry1 entry2 ...]

Built for the orbit_transfer cleanup survey (2026-07-26): produce the execution
paths, hubs, entry points and prune candidates for a campaign without reading
every file. Text-based, so it is fast and needs no MATLAB licence.

METHOD AND ITS LIMITS -- read before trusting the output:
  * Comments are stripped with a crude unquoted-'%' scan. Good enough for this
    codebase's style; it is not a MATLAB parser.
  * A call is "name appears as a word". That OVER-counts (a name mentioned in a
    string, or a variable that shadows a function name) and cannot see dynamic
    dispatch (feval, str2func, function handles built from strings).
  * It scans ONE folder tree. Cross-campaign callers are INVISIBLE, so the
    "orphan" list is a list of CANDIDATES, never a verdict. Always confirm with
    a repo-wide `grep -rlw <name> --include='*.m' .` before deleting anything.
    In the earth survey this mattered: hamiltonian_const_check looked orphaned
    here and is in fact called by the CR3BP campaign.
"""
import os, re, sys, json, hashlib, collections


def strip_comments(src):
    out = []
    for line in src.split('\n'):
        quoted, cut = False, len(line)
        for k, ch in enumerate(line):
            if ch == "'":
                quoted = not quoted
            elif ch == '%' and not quoted:
                cut = k
                break
        out.append(line[:cut])
    return '\n'.join(out)


def build(root):
    # NOTE: keyed by basename, so DUPLICATE basenames collapse to one entry and
    # their callers are merged onto whichever copy os.walk saw last. That
    # silently misattributes hubs in any tree that vendors copies -- it did
    # exactly that on GTO_tulip, where PSR/lib duplicates 14 files verbatim.
    # duplicates() below reports them; always read that section first.
    files, dupes = {}, collections.defaultdict(list)
    for dp, _, fns in os.walk(root):
        for fn in fns:
            if fn.endswith('.m'):
                dupes[fn[:-2]].append(os.path.join(dp, fn))
                files[fn[:-2]] = os.path.join(dp, fn)
    build.dupes = {n: ps for n, ps in dupes.items() if len(ps) > 1}
    names = set(files)
    calls = collections.defaultdict(set)
    for n, p in files.items():
        src = strip_comments(open(p, errors='ignore').read())
        for other in names:
            if other == n:
                continue
            if re.search(r'\b' + re.escape(other) + r'\b', src):
                calls[n].add(other)
    callers = collections.defaultdict(set)
    for a, bs in calls.items():
        for b in bs:
            callers[b].add(a)
    area = {n: (os.path.dirname(files[n]).replace(root, '').strip('/') or 'root') for n in names}
    return files, names, calls, callers, area


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    root = sys.argv[1].rstrip('/')
    entries = sys.argv[2:]
    files, names, calls, callers, area = build(root)
    tests = {n for n in names if n.startswith('test_')}
    nontest = names - tests
    UTIL = {'module_root', 'optdef', 'setup_paths'}

    entry = [n for n in sorted(nontest) if not (callers[n] - tests)]
    orphan = [n for n in entry if not calls[n]]
    print('=== %s: %d non-test, %d test ===' % (root, len(nontest), len(tests)))
    if build.dupes:
        print('\n-- !! DUPLICATE BASENAMES (%d): hub counts below are MERGED across'
              ' these copies and are not trustworthy for them --' % len(build.dupes))
        for n, ps in sorted(build.dupes.items()):
            hs = {hashlib.md5(open(q, 'rb').read()).hexdigest() for q in ps}
            print('   %-26s %-10s %s' % (n, 'IDENTICAL' if len(hs) == 1 else 'DIFFER',
                                         ' | '.join(q.replace(root + '/', '') for q in ps)))
    print('\n-- ENTRY POINTS (%d) --' % len(entry))
    for n in entry:
        print('   %-34s %-20s calls %d' % (n, area[n], len(calls[n])))
    print('\n-- ORPHAN CANDIDATES (%d) -- confirm repo-wide before deleting --' % len(orphan))
    for n in orphan:
        print('   %-34s %s' % (n, area[n]))
    print('\n-- HUBS --')
    for c, n in sorted(((len(callers[n] - tests), n) for n in nontest), reverse=True)[:12]:
        print('   %-34s %-20s called by %d' % (n, area[n], c))
    untested = [n for n in nontest if not (callers[n] & tests)]
    print('\n-- NOT REACHED BY ANY TEST: %d of %d --' % (len(untested), len(nontest)))

    def tree(rt, depth=0, seen=None, lines=None, maxd=3):
        if seen is None:
            seen, lines = set(), []
        lines.append('%s%s  [%s]' % ('  ' * depth, rt, area.get(rt, '?')))
        if depth >= maxd or rt in seen:
            return lines
        seen = seen | {rt}
        for k in sorted(c for c in calls.get(rt, []) if c not in UTIL and not c.startswith('test_')):
            tree(k, depth + 1, seen, lines, maxd)
        return lines

    for e in entries:
        print('\n' + '=' * 62 + '\nEXECUTION PATH: %s\n' % e + '=' * 62)
        print('\n'.join(tree(e)))

    out = '/tmp/callgraph_%s.json' % os.path.basename(root)
    json.dump({'calls': {k: sorted(v) for k, v in calls.items()},
               'callers': {k: sorted(v) for k, v in callers.items()},
               'area': area}, open(out, 'w'))
    print('\ngraph -> %s' % out)


if __name__ == '__main__':
    main()
