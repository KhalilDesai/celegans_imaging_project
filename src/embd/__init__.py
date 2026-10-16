"""Python 3 port of the EMBD downstream analysis (Green et al. 2024, Cell).

Faithful port: module and function names follow the original files in development/ and
embryos_analyze/, so every ported function can be diffed against its original. Behavioural
differences are listed in docs/NUMERICAL_DIFFERENCES.md and docs/DECISIONS_AND_ISSUES.md.

The input is the per-gene parameter table (docs/PARAMETER_TABLE.md). The image-analysis code
that produces it is not ported yet; calling those paths raises
NotImplementedError("image analysis not yet ported").
"""
__version__ = '0.1.0.dev0'
