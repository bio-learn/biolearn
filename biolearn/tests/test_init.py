import biolearn


def test_version_number():
    try:
        assert biolearn.__version__ == biolearn._version.__version__
    except AttributeError:
        assert biolearn.__version__ == "0+unknown"
