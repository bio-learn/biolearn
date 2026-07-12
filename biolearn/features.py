"""Stable description of the features a model requires, plus the validator
that enforces them.

Introduced in the 1.0 line. The return type of ``Model.required_features()``
and ``MissingFeaturesError`` are stable public API: breaking either requires a
major version bump.
"""

from collections.abc import Mapping
from dataclasses import dataclass

VALID_LAYERS = (
    "dnam",
    "rna",
    "protein_olink",
    "protein_alamar",
    "clinical",
)


@dataclass(frozen=True)
class RequiredFeatures(Mapping):
    """What a clock consumes, as a frozen value.

    Behaves as the documented mapping ``{"layer", "features", "metadata"}``
    (so ``rf["features"]`` and ``dict(rf)`` work) and as an object with
    ``.layer`` / ``.features`` / ``.metadata`` attributes. All three keys are
    always present; ``features`` and ``metadata`` default to empty tuples.
    """

    layer: str
    features: tuple = ()
    metadata: tuple = ()

    def __getitem__(self, key):
        return {
            "layer": self.layer,
            "features": list(self.features),
            "metadata": list(self.metadata),
        }[key]

    def __iter__(self):
        return iter(("layer", "features", "metadata"))

    def __len__(self):
        return 3


class MissingFeaturesError(ValueError):
    """Raised when a GeoData lacks features or metadata a clock requires.

    Subclasses ``ValueError`` so existing ``except ValueError`` handlers keep
    working while callers can catch this specific case.
    """


# layer -> (GeoData attribute, axis carrying feature names). The ONLY place
# that knows the index-vs-columns asymmetry: dnam/rna/protein are
# features-as-rows (features in .index); clinical is samples-as-rows
# (biomarkers in .columns), matching metadata.
_LAYER_FRAME = {
    "dnam": ("dnam", "index"),
    "rna": ("rna", "index"),
    "protein_olink": ("protein_olink", "index"),
    "protein_alamar": ("protein_alamar", "index"),
    "clinical": ("clinical", "columns"),
}

_LAYER_NOUN = {
    "dnam": "CpG sites",
    "rna": "genes",
    "protein_olink": "proteins",
    "protein_alamar": "proteins",
    "clinical": "clinical markers",
}

_MISSING_PREVIEW_LIMIT = 5


def _available_features(geo_data, layer):
    attr, axis = _LAYER_FRAME[layer]
    frame = getattr(geo_data, attr, None)
    if frame is None:
        return set()
    return set(getattr(frame, axis))


def validate_required_features(model, geo_data):
    """Raise MissingFeaturesError if geo_data lacks what the model requires."""
    req = model.required_features()

    available = _available_features(geo_data, req.layer)
    missing_features = [f for f in req.features if f not in available]

    metadata = getattr(geo_data, "metadata", None)
    metadata_cols = set(metadata.columns) if metadata is not None else set()
    missing_metadata = [m for m in req.metadata if m not in metadata_cols]

    if missing_features or missing_metadata:
        raise MissingFeaturesError(
            _format_missing(model, req, missing_features, missing_metadata)
        )


def _format_missing(model, req, missing_features, missing_metadata):
    details = getattr(model, "details", None) or {}
    name = details.get("name") if isinstance(details, dict) else None
    label = f" for model '{name}'" if name else ""
    messages = []

    if missing_features:
        missing = sorted(missing_features)
        noun = _LAYER_NOUN.get(req.layer, "features")
        preview = ", ".join(missing[:_MISSING_PREVIEW_LIMIT])
        remaining = len(missing) - _MISSING_PREVIEW_LIMIT
        if remaining > 0:
            preview = (
                f"showing first {_MISSING_PREVIEW_LIMIT}: {preview} "
                f"(+{remaining} more)"
            )
        messages.append(
            f"Missing required {noun}{label} "
            f"({len(missing)}/{len(req.features)}): {preview}. "
            f"Provide {req.layer} data with these {noun} or use an "
            f"imputation method that includes them."
        )

    if missing_metadata:
        cols = ", ".join(sorted(missing_metadata))
        messages.append(
            f"Missing required metadata{label}: {cols}. "
            f"Provide these columns in geo_data.metadata."
        )

    return " ".join(messages)
