import numpy as np
import pandas as pd
from lifelines.utils import concordance_index

from biolearn.mortality import calculate_c_index


class _StubData:
    def __init__(self, metadata):
        self.metadata = metadata


def _make_data(n=200, seed=7):
    rng = np.random.default_rng(seed)
    ids = [f"s{i}" for i in range(n)]
    age = rng.uniform(50, 90, size=n)
    # mortality driven by age plus noise
    years_until_death = np.clip(
        (100 - age) / 4 + rng.normal(0, 2, size=n), 0.1, None
    )
    dead = rng.uniform(size=n) < 0.6
    metadata = pd.DataFrame(
        {"age": age, "dead": dead.astype(int), "years_until_death": years_until_death},
        index=ids,
    )
    predictors = pd.DataFrame(
        {
            # a clock that is essentially chronological age (plus assay noise):
            # predictive of mortality, but adds nothing beyond age. Small noise
            # keeps the age-regression residuals well-defined; an EXACT copy
            # leaves only float-epsilon residuals that are still monotone in age.
            "AgeCopy": age + rng.normal(0, 0.3, size=n),
            # a clock with genuine age-independent signal
            "Informative": (100 - years_until_death * 4) + rng.normal(0, 1, size=n),
        },
        index=ids,
    )
    return _StubData(metadata), predictors


def test_default_matches_direct_concordance_and_columns():
    data, predictors = _make_data()
    result = calculate_c_index(data, predictors)
    assert list(result.columns) == ["Clock", "C_index"]
    expected = concordance_index(
        event_times=data.metadata["years_until_death"],
        predicted_scores=-predictors["AgeCopy"].astype(float),
        event_observed=data.metadata["dead"],
    )
    got = result.loc[result["Clock"] == "AgeCopy", "C_index"].iloc[0]
    assert abs(got - expected) < 1e-12


def test_bootstrap_ci_brackets_estimate_and_is_deterministic():
    data, predictors = _make_data()
    result = calculate_c_index(data, predictors, ci_bootstrap_samples=200)
    assert {"CI95_low", "CI95_high"} <= set(result.columns)
    for _, row in result.iterrows():
        assert row["CI95_low"] < row["C_index"] < row["CI95_high"]
    again = calculate_c_index(data, predictors, ci_bootstrap_samples=200)
    pd.testing.assert_frame_equal(result, again)


def test_ci_narrows_with_more_subjects():
    small_data, small_pred = _make_data(n=60, seed=3)
    large_data, large_pred = _make_data(n=600, seed=3)
    small = calculate_c_index(small_data, small_pred, ci_bootstrap_samples=200)
    large = calculate_c_index(large_data, large_pred, ci_bootstrap_samples=200)
    width = lambda df: (df["CI95_high"] - df["CI95_low"]).iloc[0]
    assert width(large) < width(small)


def test_age_adjustment_removes_age_only_signal():
    data, predictors = _make_data(n=400, seed=11)
    unadjusted = calculate_c_index(data, predictors)
    adjusted = calculate_c_index(data, predictors, adjust_for_age=True)

    def value(df, clock):
        return df.loc[df["Clock"] == clock, "C_index"].iloc[0]

    # the age-copy clock predicts mortality before adjustment...
    assert value(unadjusted, "AgeCopy") > 0.6
    # ...and collapses to chance once age is regressed out
    assert abs(value(adjusted, "AgeCopy") - 0.5) < 0.05
    # while the genuinely informative clock retains signal beyond age
    assert value(adjusted, "Informative") > 0.55
