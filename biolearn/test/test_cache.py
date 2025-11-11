import os
import shutil
import time
import pytest
import pickle

from biolearn.cache import LocalFolderCache

CATEGORY_EXPR = "expr"
CATEGORY_META = "meta"
VERSION_V1 = "v1"
VERSION_V2 = "v2"


@pytest.fixture
def cache(tmp_path):
    """Create an isolated LocalFolderCache per test."""
    cache_path = tmp_path / "cache"
    max_size_gb = 0.000000954  # 1 KB
    c = LocalFolderCache(str(cache_path), max_size_gb)
    yield c
    # Cleanup handled by tmp_path fixture


def test_store_and_get(cache):
    cache.store("key1", "value1", CATEGORY_EXPR, VERSION_V1)
    assert cache.get("key1", CATEGORY_EXPR, VERSION_V1) == "value1"


def test_get_nonexistent_key(cache):
    assert cache.get("missing", CATEGORY_EXPR, VERSION_V1) is None


def test_store_and_get_multiple_keys(cache):
    cache.store("key1", "value1", CATEGORY_EXPR, VERSION_V1)
    cache.store("key2", "value2", CATEGORY_EXPR, VERSION_V1)
    assert cache.get("key1", CATEGORY_EXPR, VERSION_V1) == "value1"
    assert cache.get("key2", CATEGORY_EXPR, VERSION_V1) == "value2"


def test_cache_size_limit(cache):
    cache.store("key1", "a" * 500, CATEGORY_EXPR, VERSION_V1)  # 500 B
    time.sleep(0.01)
    cache.store("key2", "b" * 800, CATEGORY_EXPR, VERSION_V1)  # 800 B
    time.sleep(0.01)
    cache.store("key3", "c" * 150, CATEGORY_EXPR, VERSION_V1)  # 150 B
    # key1 should be evicted (LRU)
    assert cache.get("key1", CATEGORY_EXPR, VERSION_V1) is None
    assert cache.get("key2", CATEGORY_EXPR, VERSION_V1) == "b" * 800
    assert cache.get("key3", CATEGORY_EXPR, VERSION_V1) == "c" * 150


def test_lru_respects_get_calls(cache):
    cache.store("key1", "a" * 500, CATEGORY_EXPR, VERSION_V1)
    time.sleep(0.01)
    cache.store("key2", "b" * 400, CATEGORY_EXPR, VERSION_V1)
    time.sleep(0.01)
    cache.get("key1", CATEGORY_EXPR, VERSION_V1)  # refresh key1
    time.sleep(0.01)
    cache.store("key3", "c" * 150, CATEGORY_EXPR, VERSION_V1)

    assert cache.get("key1", CATEGORY_EXPR, VERSION_V1) == "a" * 500
    assert cache.get("key2", CATEGORY_EXPR, VERSION_V1) is None
    assert cache.get("key3", CATEGORY_EXPR, VERSION_V1) == "c" * 150


def test_clear(cache):
    cache.store("key1", "value1", CATEGORY_EXPR, VERSION_V1)
    cache.store("key2", "value2", CATEGORY_EXPR, VERSION_V1)
    cache.clear()
    assert cache.get("key1", CATEGORY_EXPR, VERSION_V1) is None
    assert cache.get("key2", CATEGORY_EXPR, VERSION_V1) is None


def test_store_item_larger_than_cache_size(cache):
    cache.store("small_key", "small", CATEGORY_EXPR, VERSION_V1)
    large_item = "x" * 2000  # 2000 B
    cache.store("large_key", large_item, CATEGORY_EXPR, VERSION_V1)
    assert cache.get("large_key", CATEGORY_EXPR, VERSION_V1) is None
    assert cache.get("small_key", CATEGORY_EXPR, VERSION_V1) == "small"


def test_remove_key(cache):
    cache.store("key1", "value1", CATEGORY_EXPR, VERSION_V1)
    cache.store("key2", "value2", CATEGORY_EXPR, VERSION_V1)
    cache.store("key3", "value3", CATEGORY_EXPR, VERSION_V1)
    cache.remove("key2")
    assert cache.get("key2", CATEGORY_EXPR, VERSION_V1) is None
    assert cache.get("key1", CATEGORY_EXPR, VERSION_V1) == "value1"
    assert cache.get("key3", CATEGORY_EXPR, VERSION_V1) == "value3"


def test_version_mismatch_returns_none(cache):
    cache.store("k", "v", CATEGORY_EXPR, VERSION_V1)
    assert cache.get("k", CATEGORY_EXPR, VERSION_V2) is None


def test_category_mismatch_returns_none(cache):
    cache.store("k", "v", CATEGORY_EXPR, VERSION_V1)
    assert cache.get("k", CATEGORY_META, VERSION_V1) is None


def test_filename_composition(cache):
    cache.store("k", "v", CATEGORY_EXPR, VERSION_V1)
    files = os.listdir(cache.path)
    assert any(f.startswith("k__expr__v1") for f in files)


def test_cleanup_removes_stale_files(cache):
    # Manually create a stale file
    stale_path = os.path.join(cache.path, "k__expr__v0.pkl")
    with open(stale_path, "wb") as f:
        f.write(b"old")
    # Trigger cleanup via normal store
    cache.store("fresh", "v", CATEGORY_EXPR, VERSION_V1)
    assert "k__expr__v0.pkl" not in os.listdir(cache.path)


def test_lru_still_works_with_versions(cache):
    cache.store("k1", "a" * 500, CATEGORY_EXPR, VERSION_V1)
    time.sleep(0.01)
    cache.store("k2", "b" * 400, CATEGORY_EXPR, VERSION_V1)
    time.sleep(0.01)
    cache.get("k1", CATEGORY_EXPR, VERSION_V1)
    time.sleep(0.01)
    cache.store("k3", "c" * 150, CATEGORY_EXPR, VERSION_V1)
    assert cache.get("k1", CATEGORY_EXPR, VERSION_V1) == "a" * 500
    assert cache.get("k2", CATEGORY_EXPR, VERSION_V1) is None
    assert cache.get("k3", CATEGORY_EXPR, VERSION_V1) == "c" * 150


def test_mixed_categories_independent_versions(cache):
    cache.store("k1", "expr_v1", CATEGORY_EXPR, VERSION_V1)
    cache.store("k1", "meta_v1", CATEGORY_META, VERSION_V1)
    # now bump expr version to v2
    assert cache.get("k1", CATEGORY_EXPR, VERSION_V2) is None
    assert cache.get("k1", CATEGORY_META, VERSION_V1) == "meta_v1"


def test_no_default_version_argument(cache):
    with pytest.raises(TypeError):
        cache.get("x", CATEGORY_EXPR)  # missing version arg


def test_legacy_files_removed_on_init(tmp_path):
    """Ensure legacy <key>.pkl files are deleted when the cache is instantiated."""
    legacy_dir = tmp_path / "cache"
    legacy_dir.mkdir(parents=True, exist_ok=True)

    # Create a pre-versioning legacy file
    legacy_file = legacy_dir / "oldkey.pkl"
    legacy_file.write_bytes(b"obsolete")

    # Confirm it exists before starting
    assert legacy_file.exists()

    # Instantiate the cache (runs _remove_legacy_files)
    LocalFolderCache(str(legacy_dir), 0.000000954)  # 1 KB limit

    # Legacy file should have been removed automatically
    assert not legacy_file.exists()


def test_versioned_files_not_removed_on_init(tmp_path):
    """Correctly versioned files survive the one-time cleanup."""
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()

    # Write a properly versioned file
    data = "keep me"
    blob = pickle.dumps(data)
    versioned = cache_dir / "goodkey__expr__v1.pkl"
    versioned.write_bytes(blob)
    assert versioned.exists()

    # Init cache (should *not* delete this file)
    LocalFolderCache(str(cache_dir), max_size_gb=1e-6)
    assert versioned.exists()
