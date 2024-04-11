import os
import shutil
import pickle

class NoCache:
    def get(self, key):
        return None

    def store(self, key, value):
        pass

    def clear(self):
        pass

class LocalFolderCache:
    def __init__(self, path, max_size_gb):
        self.path = path
        self.max_size_bytes = max_size_gb * 1024 * 1024 * 1024  # Convert GB to bytes
        os.makedirs(self.path, exist_ok=True)

    def get(self, key):
        file_path = os.path.join(self.path, key)
        if os.path.exists(file_path):
            try:
                with open(file_path, 'rb') as file:
                    value = pickle.load(file)
                return value
            except (pickle.UnpicklingError, EOFError, FileNotFoundError, IOError):
                return None
        return None

    def store(self, key, value):
        file_path = os.path.join(self.path, key)
        serialized_value = pickle.dumps(value)
        value_size = len(serialized_value)
        if value_size > self.max_size_bytes:
            return
        with open(file_path, 'wb') as file:
            file.write(serialized_value)
        self._cleanup()

    def _cleanup(self):
        total_size = sum(os.path.getsize(os.path.join(self.path, f)) for f in os.listdir(self.path))
        if total_size > self.max_size_bytes:
            files = sorted(os.listdir(self.path), key=lambda f: os.path.getatime(os.path.join(self.path, f)), reverse=True)
            while total_size > self.max_size_bytes and files:
                file_to_remove = files.pop()
                file_path = os.path.join(self.path, file_to_remove)
                total_size -= os.path.getsize(file_path)
                os.remove(file_path)

    def clear(self):
        # Remove the cache directory and all its contents
        shutil.rmtree(self.path, ignore_errors=True)