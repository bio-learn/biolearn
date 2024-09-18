import os
import shutil
import pickle


class NoCache:
    """
    A cache implementation that does nothing.
    """

    def get(self, key):
        """
        Retrieves an item from the cache.

        Args:
            key (str): The key for which the cache value is retrieved.

        Returns:
            None: Always returns None, simulating a cache miss.
        """
        return None

    def store(self, key, value):
        """
        Stores an item in the cache.

        Args:
            key (str): The key under which the value is stored.
            value (any): The value to be stored.
        """
        pass

    def clear(self):
        """
        Clears all items from the cache.
        """
        pass

    def remove(self, key):
        """
        Remove a specific entry from the cache
        """
        pass


class LocalFolderCache:
    """
    A cache that stores data as serialized objects in a local directory. Implements basic LRU (Least Recently Used)
    cleanup to maintain the cache within a specified maximum size limit.
    """

    def __init__(self, path, max_size_gb):
        """
        Initializes the LocalFolderCache instance.

        Args:
            path (str): The filesystem path to the directory to be used for cache storage.
            max_size_gb (float): The maximum size of the cache in gigabytes.
        """
        self.path = path
        self.max_size_bytes = (
            max_size_gb * 1024 * 1024 * 1024
        )  # Convert GB to bytes
        os.makedirs(self.path, exist_ok=True)

    def get(self, key):
        """
        Retrieves an item from the cache by its key, if it exists.

        Args:
            key (str): The key for which the cache value is retrieved.

        Returns:
            any: The value stored in the cache, or None if the key does not exist or an error occurs.
        """
        file_path = os.path.join(self.path, key)
        if os.path.exists(file_path):
            try:
                with open(file_path, "rb") as file:
                    value = pickle.load(file)
                return value
            except (
                pickle.UnpicklingError,
                EOFError,
                FileNotFoundError,
                IOError,
            ):
                return None
        return None

    def store(self, key, value):
        """
        Stores an item in the cache under the specified key, managing the cache size to stay within limits.

        Args:
            key (str): The key under which the value is stored.
            value (any): The value to be stored.
        """
        file_path = os.path.join(self.path, key)
        serialized_value = pickle.dumps(value)
        value_size = len(serialized_value)
        if value_size > self.max_size_bytes:
            return
        with open(file_path, "wb") as file:
            file.write(serialized_value)
        self._cleanup()

    def _cleanup(self):
        """
        Cleans up the cache directory to maintain it within the maximum size limit by removing the least
        recently accessed files first. Called automatically everytime a new item is stored.
        """
        total_size = sum(
            os.path.getsize(os.path.join(self.path, f))
            for f in os.listdir(self.path)
        )
        if total_size > self.max_size_bytes:
            files = sorted(
                os.listdir(self.path),
                key=lambda f: os.path.getatime(os.path.join(self.path, f)),
                reverse=True,
            )
            while total_size > self.max_size_bytes and files:
                file_to_remove = files.pop()
                file_path = os.path.join(self.path, file_to_remove)
                total_size -= os.path.getsize(file_path)
                os.remove(file_path)

    def clear(self):
        """
        Clears the cache by removing the cache directory and all its contents.
        """
        shutil.rmtree(self.path, ignore_errors=True)

    def remove(self, key):
        """
        Remove a specific entry from the cache
        """
        file_path = os.path.join(self.path, key)
        os.remove(file_path)
