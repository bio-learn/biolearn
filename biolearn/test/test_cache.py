import os
import shutil
import pytest
import time
from biolearn.cache import LocalFolderCache


@pytest.fixture
def cache():
    cache_path = "test_cache"
    max_size_gb = 0.000000954  # 1 KB for testing purposes
    cache = LocalFolderCache(cache_path, max_size_gb)
    yield cache
    # Clean up the cache directory after each test
    shutil.rmtree(cache_path, ignore_errors=True)


def test_store_and_get(cache):
    cache.store("key1", "value1")
    assert cache.get("key1") == "value1"


def test_get_nonexistent_key(cache):
    assert cache.get("nonexistent_key") is None


def test_store_and_get_multiple_keys(cache):
    cache.store("key1", "value1")
    cache.store("key2", "value2")
    assert cache.get("key1") == "value1"
    assert cache.get("key2") == "value2"


def test_cache_size_limit(cache):
    # Store items exceeding the cache size limit
    cache.store("key1", "a" * 500)  # 500 bytes
    time.sleep(0.1)  # Ensure file timestamps are different
    cache.store("key2", "b" * 800)  # 800 bytes
    time.sleep(0.1)
    cache.store("key3", "c" * 150)  # 150 bytes

    # Check that the least recently used item (key1) is evicted
    assert cache.get("key1") is None
    assert cache.get("key2") == "b" * 800
    assert cache.get("key3") == "c" * 150


def test_lru_respects_get_calls(cache):
    # Store items exceeding the cache size limit
    cache.store("key1", "a" * 500)
    time.sleep(0.1)  # Ensure file timestamps are different
    cache.store("key2", "b" * 400)
    time.sleep(0.1)
    cache.get("key1")  # Ensure Key1 most recently used
    time.sleep(0.1)
    cache.store("key3", "c" * 150)

    # Check that the least recently used item (key1) is evicted
    assert cache.get("key1") == "a" * 500
    assert cache.get("key2") is None
    assert cache.get("key3") == "c" * 150


def test_clear(cache):
    cache.store("key1", "value1")
    cache.store("key2", "value2")
    cache.clear()
    assert cache.get("key1") is None
    assert cache.get("key2") is None


def test_store_item_larger_than_cache_size(cache):
    # Store a small item in the cache
    small_item = "small"
    cache.store("small_key", small_item)
    time.sleep(0.1)  # Ensure file timestamps are different

    large_item = "x" * 2000  # 2000 bytes
    cache.store("large_key", large_item)

    assert cache.get("large_key") is None
    assert cache.get("small_key") == small_item


def test_remove_key(cache):
    # Store multiple keys in the cache
    cache.store("key1", "value1")
    cache.store("key2", "value2")
    cache.store("key3", "value3")

    # Remove one key
    cache.remove("key2")

    # Verify that the removed key is no longer in the cache
    assert cache.get("key2") is None

    # Verify that the other keys are still in the cache
    assert cache.get("key1") == "value1"
    assert cache.get("key3") == "value3"
