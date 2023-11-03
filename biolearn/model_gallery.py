import pandas as pd
from biolearn.model import (
    model_definitions,
    LinearMethylationModel,
    ImputationDecorator,
)
from biolearn.imputation import (
    biolearn_impute,
    hybrid_impute,
    impute_from_average,
)
from biolearn.util import get_data_file


class ModelGallery:
    """
    ModelGallery manages a collection of Models that can be run on biological data to produce predictions.
    It supports retrieving models by name with various imputation methods, and searching models based on species or tissue.
    """

    def __init__(self, models=model_definitions):
        """
        Initializes the ModelGallery instance.

        Args:
            models (dict): A dictionary of model definitions.
        """
        self.models = {}
        for name, model_def in models.items():
            if not isinstance(model_def, dict):
                raise ValueError(
                    f"Expected dictionary for model definition, got {type(model_def)} for model {name}"
                )
            self.models[name] = LinearMethylationModel.from_definition(
                model_def
            )

    def get(self, name, imputation_method="default"):
        """
        Retrieves a model by its name with specified imputation method. Possible imputation values are:
        'none', 'biolearn', 'averaging', 'dunedin', 'sesame_450k'. If no value is passed it will default to sesame_450k.

        Args:
            name (str): The name of the model.
            imputation_method (str, optional): The imputation method to use

        Returns:
            object: The requested model instance, possibly wrapped in an ImputationDecorator.

        Raises:
            KeyError: If the model with the specified name is not found.
            ValueError: If an invalid imputation method is specified.
        """
        if name not in self.models:
            raise KeyError(f"Model not found: {name}")

        # Selecting imputation method
        model_instance = self.models[name]
        if imputation_method == "none":
            return model_instance
        elif imputation_method == "biolearn":
            return ImputationDecorator(
                model_instance, lambda dnam, cpgs: biolearn_impute(dnam)
            )
        elif imputation_method == "averaging":
            return ImputationDecorator(
                model_instance,
                lambda dnam, cpgs: impute_from_average(dnam, cpgs),
            )
        elif imputation_method == "dunedin":
            # DunedinPACE method
            gold_means_file = get_data_file("DunedinPACE_Gold_Means.csv")
            df = pd.read_csv(gold_means_file, index_col=0)
            gold_averages = df["mean"]
            return ImputationDecorator(
                model_instance,
                lambda dnam, cpgs: hybrid_impute(dnam, gold_averages, cpgs),
            )
        elif imputation_method in ("sesame_450k", "default"):
            # Sesame method for 450k arrays
            sesame_450k_file = get_data_file("sesame_450k_median.csv")
            df = pd.read_csv(sesame_450k_file, index_col=0)
            gold_medians = df["median"]
            return ImputationDecorator(
                model_instance,
                lambda dnam, cpgs: hybrid_impute(dnam, gold_medians, cpgs),
            )
        else:
            raise ValueError(f"Invalid imputation method: {imputation_method}")

    def search(self, species=None, tissue=None):
        """
        Search models based on species and tissue criteria.

        Args:
            species (str, optional): The species to filter models by.

            tissue (str, optional): The tissue to filter models by.

        Returns:
            dict: A dictionary of models that match the specified criteria.
        """
        matches = {}
        for name, model in self.models.items():
            if (species is None or model.metadata["species"] == species) and (
                tissue is None or model.metadata["tissue"] == tissue
            ):
                matches[name] = model
        return matches
