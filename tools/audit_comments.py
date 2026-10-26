"""Audit: are the original's explanatory comments preserved in the port?

For every function present in both an original module and its port, collect the original's
explanatory comments (``# ...`` text and docstring lines; commented-out code is excluded by a
syntax check) and report those whose text is missing from the ported function.

Run: python tools/audit_comments.py [--verbose]
"""
import ast
import os
import re
import sys
import tokenize
import io

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAIRS = [('embryos_analyze/varLookup.py', 'src/embd/varLookup.py'),
         ('development/fileLookup.py', 'src/embd/fileLookup.py'),
         ('development/params_use_assign.py', 'src/embd/params_use_assign.py'),
         ('embryos_analyze/myFunc.py', 'src/embd/myFunc.py'),
         ('embryos_analyze/myMath.py', 'src/embd/myMath.py'),
         ('embryos_analyze/myFigure.py', 'src/embd/myFigure.py'),
         ('development/embdFunc.py', 'src/embd/embdFunc.py'),
         ('embryos_analyze/db_utils_embryos.py', 'src/embd/db_utils_embryos.py'),
         ('development/db_utils_genes.py', 'src/embd/db_utils_genes.py'),
         ('development/RNAiClass.py', 'src/embd/RNAiClass.py'),
         ('development/phenoDistCalc.py', 'src/embd/phenoDistCalc.py'),
         ('development/PADassist.py', 'src/embd/PADassist.py'),
         ('development/groupAssesment.py', 'src/embd/groupAssesment.py'),
         ('development/makeUMAPs_V2.py', 'src/embd/makeUMAPs_V2.py')]


def norm(s):
    return re.sub(r'\s+', ' ', s).strip().lower()


def is_code(text):
    t = text.strip()
    if not t:
        return True
    for candidate in (t, 'print ' + t if t.startswith('print') else t):
        try:
            ast.parse(candidate.replace('print ', 'print(', 1) + ')' if candidate.startswith('print ') else candidate)
            # bare words/sentences also parse as names; require code-ish punctuation
            if re.search(r'[=()\[\].]', t) and not re.match(r'^[A-Za-z ,\'-]+[.!?]?$', t):
                return True
        except SyntaxError:
            pass
    return False


def functions(src):
    """name -> (start_line, end_line) for top-level defs and class methods (by indentation)."""
    lines = src.splitlines()
    out, stack = {}, []
    for i, line in enumerate(lines):
        m = re.match(r'^(\s*)def (\w+)\(', line)
        if m:
            out.setdefault(m.group(2), []).append(i)
    res = {}
    for name, starts in out.items():
        s = starts[-1]
        ind = len(lines[s]) - len(lines[s].lstrip())
        e = s + 1
        while e < len(lines) and (not lines[e].strip() or len(lines[e]) - len(lines[e].lstrip()) > ind or lines[e].lstrip().startswith('#')):
            e += 1
        res[name] = '\n'.join(lines[s:e])
    return res


def comments(body):
    found = []
    for line in body.splitlines():
        if '#' in line:
            # crude but adequate: text after the first # not inside quotes
            q, cut = None, None
            for k, ch in enumerate(line):
                if ch in '\'"':
                    q = None if q == ch else (ch if q is None else q)
                elif ch == '#' and q is None:
                    cut = k
                    break
            if cut is not None:
                text = line[cut + 1:].strip()
                if text and not is_code(text):
                    found.append(text)
    m = re.search(r'(\'\'\'|""")(.*?)(\1)', body, re.S)
    if m:
        for line in m.group(2).splitlines():
            if line.strip():
                found.append(line.strip())
    return found


def main(verbose):
    total, missing_total = 0, 0
    for orig, port in PAIRS:
        o = functions(open(os.path.join(REPO, orig), encoding='latin-1').read())
        p = functions(open(os.path.join(REPO, port), encoding='utf-8').read())
        for name in sorted(set(o) & set(p)):
            pbody = norm(p[name])
            cs = comments(o[name])
            miss = [c for c in cs if norm(c)[:60] not in pbody]
            total += len(cs)
            missing_total += len(miss)
            if miss:
                print('%s :: %s  missing %d/%d' % (port, name, len(miss), len(cs)))
                if verbose:
                    for c in miss:
                        print('      - ' + c[:110])
    print('TOTAL explanatory comment lines: %d, missing in port: %d' % (total, missing_total))


if __name__ == '__main__':
    main('--verbose' in sys.argv)
