import numpy as np
import pandas as pd

from biolearn.dunedin_pace import quantile_normalize_using_target


def _target():
    return list(np.linspace(0.0, 1.0, 10))


def test_quantile_normalize_accepts_read_only_array():
    # Arrays from DataFrame.values can be read-only (e.g. pandas
    # copy-on-write in pandas >= 3.0). Normalization must not fail with
    # "assignment destination is read-only". See issue #195.
    data = pd.DataFrame(np.random.rand(10, 3)).values
    data.flags.writeable = False

    result = quantile_normalize_using_target(data, _target())

    assert result.shape == (10, 3)


def test_quantile_normalize_does_not_mutate_input():
    data = pd.DataFrame(np.random.rand(10, 3)).values
    original = data.copy()

    quantile_normalize_using_target(data, _target())

    np.testing.assert_array_equal(data, original)


def test_quantile_normalize_maps_values_onto_target_distribution():
    target = _target()
    data = pd.DataFrame(np.random.rand(10, 3)).values

    result = quantile_normalize_using_target(data, target)

    sorted_target = np.sort(target)
    midpoints = 0.5 * (sorted_target[:-1] + sorted_target[1:])
    allowed = np.round(np.concatenate([sorted_target, midpoints]), 6)
    assert np.all(np.isin(np.round(result, 6), allowed))
