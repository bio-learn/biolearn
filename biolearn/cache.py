from __future__ import annotations

import os
import pickle
import shutil
from typing import Any

__all__ = ["NoCache", "LocalFolderCache"]


def _compose_filename(key: str, category: str, version: str) -> str:
    """
    Construct the filename for a cache entry.

    Args:
        key: Unique key identifying the cache entry.
        category: Cache category name.
        version: Version string for invalidation.

    Returns:
        Filename string including key, category, and version.
    """
    return f"{key}__{category}__{version}.pkl"


class NoCache:
    """A no-op cache implementation (always misses)."""

    def get(self, key: str, category: str, version: str) -> Any:
        """
        Retrieve a cache entry.

        Args:
            key: Unique key identifying the cache entry.
            category: Cache category name.
            version: Version string for invalidation.

        Returns:
            None always, indicating a cache miss.
        """
        return None

    def store(self, key: str, value: Any, category: str, version: str) -> None:
        """
        No-op store method.

        Args:
            key: Unique key identifying the cache entry.
            value: Object to store.
            category: Cache category name.
            version: Version string for invalidation.
        """
        pass

    def clear(self) -> None:
        """
        Remove all cache entries (no-op).

        Args:
            None
        """
        pass

    def remove(self, key: str) -> None:
        """
        Remove cache entries matching the given key (no-op).

        Args:
            key: Unique key identifying entries to remove.
        """
        pass


class LocalFolderCache:
    """LRU-style cache that stores pickled objects in a local directory."""

    def __init__(self, path: str, max_size_gb: float):
        """
        Initialize the LocalFolderCache.

        Args:
            path: Directory path for storing cache files.
            max_size_gb: Maximum total cache size in gigabytes.
        """
        self.path = path
        self.max_size_bytes = int(max_size_gb * 1024**3)
        os.makedirs(self.path, exist_ok=True)
        # One‑time clean‑up: remove *legacy* cache files (no category/version).
        # This was introduced May 2025, should be reasonable removed in a year or two.
        self._remove_legacy_files()

    def get(self, key: str, category: str, version: str) -> Any | None:
        """
        Retrieve a cached value.

        Args:
            key: Unique key identifying the cache entry.
            category: Cache category name.
            version: Version string for invalidation.

        Returns:
            The cached object, or None if not found or on load error.
        """
        filepath = os.path.join(
            self.path, _compose_filename(key, category, version)
        )
        if not os.path.exists(filepath):
            return None
        try:
            with open(filepath, "rb") as f:
                return pickle.load(f)
        except (pickle.UnpicklingError, EOFError, OSError):
            return None

    def store(self, key: str, value: Any, category: str, version: str) -> None:
        """
        Store a value in the cache under the given key/category/version.

        Args:
            key: Unique key identifying the cache entry.
            value: Object to store.
            category: Cache category name.
            version: Version string for invalidation.

        Notes:
            If the serialized object exceeds the cache size limit, it will not be stored.
            Before storing, stale files in the same category but different version are removed.
        """
        data = pickle.dumps(value)
        if len(data) > self.max_size_bytes:
            return

        self._remove_stale(category, version)

        filepath = os.path.join(
            self.path, _compose_filename(key, category, version)
        )
        with open(filepath, "wb") as f:
            f.write(data)

        self._enforce_size_limit()

    def clear(self) -> None:
        """
        Remove all cache files.

        Args:
            None
        """
        shutil.rmtree(self.path, ignore_errors=True)
        os.makedirs(self.path, exist_ok=True)

    def remove(self, key: str) -> None:
        """
        Remove cache files matching the given key across categories and versions.

        Args:
            key: Unique key identifying entries to remove.
        """
        prefix = f"{key}__"
        for fname in os.listdir(self.path):
            if fname.startswith(prefix):
                os.remove(os.path.join(self.path, fname))

    def _remove_stale(self, category: str, allowed_version: str) -> None:
        """
        Delete files in the given category that do not match the allowed version.

        Args:
            category: Cache category name.
            allowed_version: Version string to retain.
        """
        for fname in os.listdir(self.path):
            parts = fname.split("__")
            if len(parts) != 3 or not fname.endswith(".pkl"):
                continue
            _, cat, ver_ext = parts
            version = ver_ext.rsplit(".", 1)[0]
            if cat == category and version != allowed_version:
                os.remove(os.path.join(self.path, fname))

    def _enforce_size_limit(self) -> None:
        """
        Evict least-recently-used files until the cache size is under the configured limit.

        Args:
            None
        """
        files = os.listdir(self.path)
        total = sum(os.path.getsize(os.path.join(self.path, f)) for f in files)
        if total <= self.max_size_bytes:
            return

        files.sort(key=lambda f: os.path.getatime(os.path.join(self.path, f)))
        while total > self.max_size_bytes and files:
            oldest = files.pop(0)
            path = os.path.join(self.path, oldest)
            total -= os.path.getsize(path)
            os.remove(path)

    def _remove_legacy_files(self) -> None:
        """Delete any file that does *not* follow the three‑part naming scheme."""
        for fname in os.listdir(self.path):
            if fname.count("__") != 2 or not fname.endswith(".pkl"):
                continue  # not a cache file
            # legacy files have <key>.pkl (zero separators) or <key>__misc.pkl (one separator)
        for fname in list(os.listdir(self.path)):
            if not fname.endswith(".pkl"):
                continue
            if fname.count("__") != 2:
                os.remove(os.path.join(self.path, fname))
