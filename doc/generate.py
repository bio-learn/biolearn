"""Generate model documentation tables for biolearn.

This script is called by the doc build to produce rst/csv fragments that
list each model's required CpG sites, coefficients, and metadata.

Patch for issue #200: models without a ``model_file`` entry (e.g.
HurdleInflammAge) now have their required CpG site count surfaced by
instantiating the model and calling ``methylation_sites()`` rather than
showing N/A.
"""

from __future__ import annotations

import importlib
import inspect
import traceback
from pathlib import Path
from typing import Optional


def get_model_cpg_sites(model_name: str, model_def: dict) -> str:
    """Return a human-readable string describing the required CpG sites.

    For models that list a ``model_file``, we report 'See coefficients file'.
    For models without one (e.g. HurdleInflammAge), we attempt to instantiate
    the model class and call ``methylation_sites()`` to count the sites
    dynamically.

    Returns a plain string suitable for inclusion in a doc table cell.
    """
    if model_def.get("model_file"):
        return "See coefficients file"

    # Try to surface sites from a live model instance.
    # We try common import paths used by biolearn models.
    candidate_modules = [
        f"biolearn.model.{model_name.lower()}",
        "biolearn.model.clocks",
        "biolearn.model.hurdle",
        "biolearn.model",
    ]

    for module_path in candidate_modules:
        try:
            mod = importlib.import_module(module_path)
            cls = getattr(mod, model_name, None)
            if cls is None or not inspect.isclass(cls):
                continue
            instance = cls()
            if not hasattr(instance, "methylation_sites"):
                continue
            sites = instance.methylation_sites()
            count = len(sites) if sites is not None else 0
            return f"{count} sites"
        except Exception:
            continue  # try next candidate

    return "N/A"


def get_model_info(model_name: str, model_def: dict) -> dict:
    """Collect documentation fields for a single model."""
    return {
        "name": model_name,
        "species": model_def.get("species", "N/A"),
        "tissue": model_def.get("tissue", "N/A"),
        "cpg_sites": get_model_cpg_sites(model_name, model_def),
        "reference": model_def.get("reference", "N/A"),
    }


def generate_model_table(models: dict) -> str:
    """Return an RST table string for all models."""
    rows = [get_model_info(name, defn) for name, defn in models.items()]

    headers = ["Model", "Species", "Tissue", "CpG Sites", "Reference"]
    col_widths = [max(len(h), max((len(str(r[k])) for r in rows), default=0))
                  for h, k in zip(headers, ["name", "species", "tissue", "cpg_sites", "reference"])]

    sep = "  ".join("-" * w for w in col_widths)
    header_row = "  ".join(h.ljust(w) for h, w in zip(headers, col_widths))

    lines = [sep, header_row, sep]
    for row in rows:
        lines.append("  ".join(
            str(row[k]).ljust(w)
            for k, w in zip(["name", "species", "tissue", "cpg_sites", "reference"], col_widths)
        ))
    lines.append(sep)
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    # Quick smoke-test: import the model registry and print the table.
    try:
        from biolearn.model.model_definitions import MODEL_DEFINITIONS
        print(generate_model_table(MODEL_DEFINITIONS))
    except ImportError as exc:
        print(f"Could not import MODEL_DEFINITIONS: {exc}")
