import pandas as pd
from biolearn.clock import clock_definitions, LinearMethylationClock, ImputationDecorator
from biolearn.imputation import biolearn_impute, hybrid_impute, impute_from_average
from biolearn.util import get_data_file


class ClockGallery:
    def __init__(self, clocks=clock_definitions):
        self.clocks = {}
        for name, clock_def in clocks.items():
            if not isinstance(clock_def, dict):
                raise ValueError(
                    f"Expected dictionary for clock definition, got {type(clock_def)} for clock {name}"
                )
            self.clocks[name] = LinearMethylationClock.from_definition(clock_def)

    def get(self, name, imputation_method="default"):
        """Get a clock by its name."""
        if name not in self.clocks:
            raise KeyError(f"Clock not found: {name}")

        # Wrap the clock in an ImputationDecorator based on the imputation_method parameter
        clock_instance = self.clocks[name]
        if imputation_method == "none":
            return clock_instance
        elif imputation_method == "biolearn":
            return ImputationDecorator(clock_instance, lambda dnam, cpgs: biolearn_impute(dnam))
        elif imputation_method == "averaging":
            return ImputationDecorator(clock_instance, lambda dnam, cpgs: impute_from_average(dnam, cpgs))
        elif imputation_method == "dunedin" or imputation_method == "default":
            gold_means_file = get_data_file("DunedinPACE_Gold_Means.csv")
            df = pd.read_csv(gold_means_file, index_col=0)
            gold_averages = df["mean"]
            return ImputationDecorator(clock_instance, lambda dnam, cpgs: hybrid_impute(dnam, gold_averages, cpgs))
        else:
            raise ValueError(f"Invalid imputation method: {imputation_method}")

    def search(self, species=None, tissue=None):
        matches = {}
        for name, clock in self.clocks.items():
            if (species is None or clock.metadata["species"] == species) and (
                tissue is None or clock.metadata["tissue"] == tissue
            ):
                matches[name] = clock
        return matches
