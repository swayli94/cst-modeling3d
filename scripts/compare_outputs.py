#!/usr/bin/env python
"""Compare example outputs of two runs, e.g. before/after a refactor.

Usage: python scripts/compare_outputs.py <baseline_dir> <new_dir>

Both directories must be example-tree copies in which
scripts/run_all_examples.sh has been executed.

Text outputs (.dat/.xyz/.igs/.txt/.in) are compared byte-wise first, then
numerically (float tokens, rel/abs tolerance). PNGs are byte-compared.
Files present only in <new_dir> are reported separately at the end.
"""
import os
import sys
import math

TOL_ABS = 1e-9
TOL_REL = 1e-7

TEXT_EXT = ('.dat', '.xyz', '.igs', '.txt', '.in')
CMP_EXT = TEXT_EXT + ('.png',)

def tokens_as_floats(path):
    floats = []
    with open(path, 'r', errors='replace') as f:
        for tok in f.read().split():
            try:
                floats.append(float(tok))
            except ValueError:
                floats.append(None)
    return floats

def compare_text(a, b):
    with open(a, 'rb') as fa, open(b, 'rb') as fb:
        if fa.read() == fb.read():
            return 'identical', 0.0
    ta, tb = tokens_as_floats(a), tokens_as_floats(b)
    if len(ta) != len(tb):
        return 'DIFFERENT (token count %d vs %d)' % (len(ta), len(tb)), math.inf
    max_diff = 0.0
    n_diff = 0
    for xa, xb in zip(ta, tb):
        if xa is None or xb is None:
            if xa != xb:
                return 'DIFFERENT (non-numeric tokens)', math.inf
            continue
        d = abs(xa - xb)
        if d > TOL_ABS and d > TOL_REL * max(abs(xa), abs(xb)):
            n_diff += 1
            max_diff = max(max_diff, d)
    if n_diff == 0:
        return 'close (formatting only)', 0.0
    return 'NUMERIC-DIFF n=%d' % n_diff, max_diff

def collect(base):
    files = {}
    for root, _, names in os.walk(base):
        for fn in names:
            if fn == '_results.txt':
                continue
            ext = os.path.splitext(fn)[1]
            if ext in CMP_EXT:
                p = os.path.join(root, fn)
                files[os.path.relpath(p, base)] = p
    return files

def main():
    base, new = sys.argv[1], sys.argv[2]
    fa, fb = collect(base), collect(new)
    n_bad = 0
    for rel in sorted(fa):
        pa = fa[rel]
        pb = fb.get(rel)
        if pb is None:
            print('%-58s MISSING in new run' % rel)
            n_bad += 1
            continue
        ext = os.path.splitext(rel)[1]
        if ext == '.png':
            with open(pa, 'rb') as f1, open(pb, 'rb') as f2:
                same = f1.read() == f2.read()
            if not same:
                print('%-58s png-bytes-differ' % rel)
                n_bad += 1
            continue
        status, md = compare_text(pa, pb)
        if not (status.startswith('identical') or status.startswith('close')):
            print('%-58s %-40s maxdiff=%.3e' % (rel, status, md))
            n_bad += 1
    only_new = sorted(set(fb) - set(fa))
    print('----')
    print('%d files compared, %d with differences, %d only in new run'
          % (len(fa), n_bad, len(only_new)))
    for rel in only_new:
        print('  only in new: %s' % rel)

if __name__ == '__main__':
    main()
