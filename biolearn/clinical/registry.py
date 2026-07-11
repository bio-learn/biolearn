"""Biomarker registry defining canonical names, units, valid ranges, and conversions.

The canonical units match NHANES conventions established in biolearn.load.
All clinical clocks expect data in these units.
"""

_REGISTRY = {
    "albumin": {
        # NHANES exposes albumin via LBDSALSI in g/L, so we use that as the
        # canonical unit. Typical human serum albumin is ~35-50 g/L.
        "unit": "g/L",
        "range": (10.0, 60.0),
        "description": "Serum albumin",
        "conversions": {
            "g/dL": lambda x: x * 10.0,
        },
    },
    "creatinine": {
        # NHANES exposes creatinine via LBDSCRSI in umol/L. Typical adult range
        # is ~50-110 umol/L.
        "unit": "umol/L",
        "range": (10.0, 1500.0),
        "description": "Serum creatinine",
        "conversions": {
            "mg/dL": lambda x: x * 88.42,
        },
    },
    "glucose": {
        "unit": "mmol/L",
        "range": (1.0, 40.0),
        "description": "Fasting glucose",
        "conversions": {
            "mg/dL": lambda x: x * 0.05551,
        },
    },
    "c_reactive_protein": {
        "unit": "mg/dL",
        "range": (0.01, 30.0),
        "description": "C-reactive protein",
        "conversions": {
            "mg/L": lambda x: x / 10.0,
            "nmol/L": lambda x: x / 95.24,
        },
    },
    "white_blood_cell_count": {
        "unit": "1000 cells/uL",
        "range": (1.0, 50.0),
        "description": "White blood cell count",
        "conversions": {},
    },
    "lymphocyte_percent": {
        "unit": "%",
        "range": (1.0, 80.0),
        "description": "Lymphocyte percentage",
        "conversions": {},
    },
    "red_blood_cell_distribution_width": {
        "unit": "%",
        "range": (8.0, 30.0),
        "description": "Red blood cell distribution width",
        "conversions": {},
    },
    "mean_cell_volume": {
        "unit": "fL",
        "range": (50.0, 130.0),
        "description": "Mean corpuscular volume",
        "conversions": {},
    },
    "alkaline_phosphate": {
        "unit": "U/L",
        "range": (10.0, 500.0),
        "description": "Alkaline phosphatase",
        "conversions": {},
    },
    "hdl_cholesterol": {
        "unit": "mmol/L",
        "range": (0.2, 5.0),
        "description": "HDL cholesterol",
        "conversions": {
            "mg/dL": lambda x: x / 38.67,
        },
    },
    "hemoglobin": {
        "unit": "g/dL",
        "range": (4.0, 22.0),
        "description": "Hemoglobin",
        "conversions": {
            "g/L": lambda x: x / 10.0,
        },
    },
    "platelet_count": {
        "unit": "1000 cells/uL",
        "range": (10.0, 1000.0),
        "description": "Platelet count",
        "conversions": {},
    },
    "mean_cell_hemoglobin": {
        "unit": "pg",
        "range": (15.0, 45.0),
        "description": "Mean corpuscular hemoglobin",
        "conversions": {},
    },
    "basophil_percent": {
        "unit": "%",
        "range": (0.0, 10.0),
        "description": "Basophil percentage",
        "conversions": {},
    },
    "lymphocyte_number": {
        "unit": "1000 cells/uL",
        "range": (0.1, 20.0),
        "description": "Lymphocyte count",
        "conversions": {},
    },
    "red_blood_cell_count": {
        "unit": "million cells/uL",
        "range": (1.0, 10.0),
        "description": "Red blood cell count",
        "conversions": {},
    },
}

# Source unit presets for common data sources.
#
# Only sources we have actually loaded and tested against are listed here. Do
# not add a preset unless you can validate it against real data from that source
# (see ``test_load_fhs_as_geodata`` for the FHS validation pattern).
_SOURCE_PRESETS = {
    # NHANES values are already in canonical units after ``load_nhanes`` runs,
    # so no conversions are needed.
    "nhanes": {},
    # FHS (frmgham2.csv) ships glucose in mg/dL; canonical is mmol/L. Other
    # FHS biomarkers (HDL/LDL/total cholesterol) aren't populated in Period 1,
    # so they're not part of the preset until we expand load_fhs to other
    # periods.
    "fhs": {
        "glucose": "mg/dL",
    },
}


class BiomarkerRegistry:
    """Registry of canonical biomarker definitions for clinical clocks.

    Provides lookup of units, valid ranges, and conversion functions
    for all biomarkers used by clinical aging clocks.
    """

    def __init__(self, registry=None):
        self._registry = registry or _REGISTRY

    def get(self, name):
        """Get biomarker definition by canonical name.

        Parameters
        ----------
        name : str
            Canonical biomarker name (e.g. 'albumin', 'creatinine').

        Returns
        -------
        dict
            Biomarker definition with keys: unit, range, description, conversions.

        Raises
        ------
        KeyError
            If the biomarker name is not in the registry.
        """
        if name not in self._registry:
            raise KeyError(
                f"Unknown biomarker: '{name}'. "
                f"Known biomarkers: {', '.join(sorted(self._registry.keys()))}"
            )
        return self._registry[name]

    def canonical_unit(self, name):
        """Return the canonical unit for a biomarker."""
        return self.get(name)["unit"]

    def valid_range(self, name):
        """Return the (min, max) valid range for a biomarker."""
        return self.get(name)["range"]

    def known_biomarkers(self):
        """Return sorted list of all known biomarker names."""
        return sorted(self._registry.keys())

    def get_source_preset(self, source_name):
        """Return unit mapping for a named data source.

        Parameters
        ----------
        source_name : str
            Source preset name (e.g. 'nhanes', 'ukbiobank').

        Returns
        -------
        dict
            Mapping of biomarker name to source unit string.
        """
        if source_name not in _SOURCE_PRESETS:
            raise ValueError(
                f"Unknown source preset: '{source_name}'. "
                f"Known presets: {', '.join(sorted(_SOURCE_PRESETS.keys()))}"
            )
        return _SOURCE_PRESETS[source_name]

    def __contains__(self, name):
        return name in self._registry

    def __len__(self):
        return len(self._registry)


BIOMARKER_REGISTRY = BiomarkerRegistry()
