"""Run-time locations (replaces the hardcoded ``FOLDER_IN = 'Z:/'``; DECISIONS, PORTING_PLAN section 7).

* ``EMBD_OUT_DIR`` stands in for the original ``Z:/`` root: inputs the original code read from
  ``Z:/`` (gene_key.csv, phenoDist CSVs) and every output live here.
* ``EMBD_SQLITE`` is the MySQL emulation database (see embd_data.mysql_emulation).
* ``EMBD_PARAMSET`` selects the parameter set returned by ``params_use_assign.get_params_use``
  (opt | gls | ms | all). This replaces editing the ``return`` line by hand (DECISIONS L5). It is an
  environment variable so that multiprocessing workers inherit it.
"""
import os


def out_dir():
    d = os.environ.get('EMBD_OUT_DIR', os.path.join(os.getcwd(), 'embd_out'))
    os.makedirs(d, exist_ok=True)
    return d


def folder_in():
    """Equivalent of ``varLookup.FOLDER_IN`` ('Z:/'), with a trailing separator."""
    return out_dir().rstrip('/') + '/'
