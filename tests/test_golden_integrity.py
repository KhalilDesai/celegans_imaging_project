"""The golden files must be exactly what MANIFEST.sha256 says (guards against accidental edits)."""
import hashlib
import os

GOLDEN = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'golden')


def test_manifest_matches_files():
    with open(os.path.join(GOLDEN, 'MANIFEST.sha256')) as f:
        entries = [line.split() for line in f if line.strip()]
    assert len(entries) == 41
    for digest, name in entries:
        with open(os.path.join(GOLDEN, name), 'rb') as g:
            assert hashlib.sha256(g.read()).hexdigest() == digest, name
