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
