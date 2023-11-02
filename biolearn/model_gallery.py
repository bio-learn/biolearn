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
        Get a model by its name.

        Args:
            name (str): The name of the model.
            imputation_method (str, optional): The imputation method to use. Defaults to "default".

        Returns:
            The requested model instance, possibly wrapped in an ImputationDecorator.

        Raises:
            KeyError: If the model with the specified name is not found.
            ValueError: If an invalid imputation method is specified.
        """
        if name not in self.models:
            raise KeyError(f"Model not found: {name}")

        # Wrap the model in an ImputationDecorator based on the imputation_method parameter
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
            gold_means_file = get_data_file("DunedinPACE_Gold_Means.csv")
            df = pd.read_csv(gold_means_file, index_col=0)
            gold_averages = df["mean"]
            return ImputationDecorator(
                model_instance,
                lambda dnam, cpgs: hybrid_impute(dnam, gold_averages, cpgs),
            )
        elif (
            imputation_method == "sesame_450k"
            or imputation_method == "default"
        ):
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
