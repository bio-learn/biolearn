from decimal import Decimal, localcontext

import numpy as np
import pandas as pd
import pytest

from biolearn.hematology import phenotypic_age

# Levine et al. (2018), Aging 10(4):573-591, Table 1 / Supplementary Table S1
# and Supplement 1 (gamma). ALP weight unrounded as in PLoS Med 16(2):e1002760.
LEVINE_COEFS = {
    "albumin": "-0.0336",
    "creatinine": "0.0095",
    "glucose": "0.1953",
    "c_reactive_protein": "0.0954",  # applied to ln(CRP in mg/dL)
    "lymphocyte_percent": "-0.0120",
    "mean_cell_volume": "0.0268",
    "red_blood_cell_distribution_width": "0.3306",
    "alkaline_phosphate": "0.00188",
    "white_blood_cell_count": "0.0554",
    "age": "0.0804",
}
LEVINE_INTERCEPT = "-19.9067"
LEVINE_GAMMA = "0.0076927"

TYPICAL = {
    "age": 35,
    "albumin": 46,
    "creatinine": 80,
    "glucose": 4.9,
    "c_reactive_protein": 0.08,
    "lymphocyte_percent": 33,
    "mean_cell_volume": 89,
    "red_blood_cell_distribution_width": 12.4,
    "alkaline_phosphate": 62,
    "white_blood_cell_count": 5.6,
}

# Mortality score so close to 1 that the published two-step formula
# saturates in float64 (1 - M is ~1e-13 here, ~1e-440 for the second row).
HIGH_RISK = [
    {
        "age": 80,
        "albumin": 30,
        "creatinine": 180,
        "glucose": 9,
        "c_reactive_protein": 2,
        "lymphocyte_percent": 12,
        "mean_cell_volume": 98,
        "red_blood_cell_distribution_width": 17,
        "alkaline_phosphate": 150,
        "white_blood_cell_count": 12,
    },
    {
        "age": 85,
        "albumin": 28,
        "creatinine": 300,
        "glucose": 15,
        "c_reactive_protein": 5,
        "lymphocyte_percent": 10,
        "mean_cell_volume": 100,
        "red_blood_cell_distribution_width": 18,
        "alkaline_phosphate": 200,
        "white_blood_cell_count": 15,
    },
]


def _published_formula(row, prec=1000):
    """Levine's two-step formula evaluated literally in high precision."""
    with localcontext() as ctx:
        ctx.prec = prec
        xb = Decimal(LEVINE_INTERCEPT)
        for name, coef in LEVINE_COEFS.items():
            value = Decimal(str(row[name]))
            if name == "c_reactive_protein":
                value = value.ln()
            xb += Decimal(coef) * value
        gamma = Decimal(LEVINE_GAMMA)
        hazard = xb.exp() * ((120 * gamma).exp() - 1) / gamma
        mortality_score = 1 - (-hazard).exp()
        return float(
            Decimal("141.50225")
            + (Decimal("-0.00553") * (1 - mortality_score).ln()).ln()
            / Decimal("0.090165")
        )


def test_issue_220_example_matches_paper_coefficients():
    # Value reported in issue #220 for the paper's coefficients.
    result = phenotypic_age(pd.DataFrame([TYPICAL]))
    assert result.iloc[0] == pytest.approx(24.38, abs=0.005)


def test_matches_published_formula_for_typical_inputs():
    rng = np.random.default_rng(0)
    n = 200
    df = pd.DataFrame(
        {
            "age": rng.uniform(20, 85, n),
            "albumin": rng.uniform(35, 52, n),
            "creatinine": rng.uniform(50, 130, n),
            "glucose": rng.uniform(4, 8, n),
            "c_reactive_protein": rng.uniform(0.02, 1.5, n),
            "lymphocyte_percent": rng.uniform(15, 45, n),
            "mean_cell_volume": rng.uniform(80, 100, n),
            "red_blood_cell_distribution_width": rng.uniform(11.5, 15, n),
            "alkaline_phosphate": rng.uniform(40, 120, n),
            "white_blood_cell_count": rng.uniform(3.5, 10, n),
        }
    )
    expected = [_published_formula(row) for _, row in df.iterrows()]

    np.testing.assert_allclose(phenotypic_age(df), expected, atol=1e-9)


def test_high_risk_inputs_are_not_capped():
    # These previously returned ~110.96 (the cap implied by clamping
    # 1 - M at 1e-5).
    result = phenotypic_age(pd.DataFrame(HIGH_RISK))
    expected = [_published_formula(row) for row in HIGH_RISK]

    assert np.all(np.isfinite(result))
    assert np.all(result > 115)
    np.testing.assert_allclose(result, expected, atol=1e-9)


def test_does_not_mutate_input():
    df = pd.DataFrame([TYPICAL] + HIGH_RISK)
    original = df.copy()

    phenotypic_age(df)

    pd.testing.assert_frame_equal(df, original)
