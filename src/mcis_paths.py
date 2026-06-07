"""Shared path resolution. Data directory comes from MCIS_DATA_DIR, falling
back to <repo>/data/ — never a machine-specific absolute path."""
import os


def repo_root():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def data_dir():
    d = os.environ.get("MCIS_DATA_DIR") or os.path.join(repo_root(), "data")
    return d.rstrip("/") + "/"


def figures_dir():
    return os.path.join(repo_root(), "figures") + "/"


def results_dir():
    return os.path.join(repo_root(), "results") + "/"
