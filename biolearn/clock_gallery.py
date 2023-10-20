import pandas as pd
from biolearn.clock import clock_definitions, LinearMethylationClock


class ClockGallery:
    def __init__(self, clocks=clock_definitions):
        self.clocks = {}
        for name, clock_def in clocks.items():
            if not isinstance(clock_def, dict):
                raise ValueError(
                    f"Expected dictionary for clock definition, got {type(clock_def)} for clock {name}"
                )
            self.clocks[name] = LinearMethylationClock.from_definition(
                clock_def
            )

    def get_by_name(self, name):
        """Get a clock by its name."""
        if name not in self.clocks:
            raise KeyError(f"Clock not found: {name}")
        return self.clocks[name]

    def search(self, species=None, tissue=None):
        matches = {}
        for name, clock in self.clocks.items():
            if (species is None or clock.metadata["species"] == species) and (
                tissue is None or clock.metadata["tissue"] == tissue
            ):
                matches[name] = clock
        return matches
