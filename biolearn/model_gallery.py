import pandas as pd
from biolearn.model import (
    model_definitions,
    LinearMethylationModel,
    LinearTranscriptomicModel,
    GrimageModel,
    SexEstimationModel,
    ImputationDecorator,
    DeconvolutionModel,
    LinearMultipartProteomicModel,
    EpiTOC2Model,
)
from biolearn.imputation import (
    hybrid_impute,
    impute_from_average,
)
from biolearn.util import get_data_file


class ModelGallery:
    """
    ModelGallery manages a collection of Models that can be run on biological data to produce predictions.
    It supports retrieving models by name with various imputation methods, and searching models based on species or tissue.
    """

    model_builders = {
        "LinearMethylationModel": LinearMethylationModel.from_definition,
        "LinearTranscriptomicModel": LinearTranscriptomicModel.from_definition,
        "GrimageModel": GrimageModel.from_definition,
        "SexEstimationModel": SexEstimationModel.from_definition,
        "DeconvolutionModel": DeconvolutionModel.from_definition,
        "LinearMultipartProteomicModel": LinearMultipartProteomicModel.from_definition,
        "EpiTOC2Model": EpiTOC2Model.from_definition,
    }

    def __init__(self, models=model_definitions):
        """
        Initializes the ModelGallery instance.

        Args:
            models (dict): A dictionary of model definitions.
        """
        self.model_definitions = {}
        for name, model_def in models.items():
            if not isinstance(model_def, dict):
                raise ValueError(
                    f"Expected dictionary for model definition, got {type(model_def)} for model {name}"
                )
            model_type = model_def["model"]["type"]
            if model_type in self.model_builders:
                self.model_definitions[name] = model_def
            else:
                raise ValueError(
                    f"Model type {model_type} does not have a known builder"
                )

    def get(self, name, imputation_method=None):
        """
        Retrieves a model by its name with a specified imputation method to handle missing data.

        Args:
            name (str): The name of the model.
            imputation_method (str, optional): The method used for imputing missing data in the model's dataset.

        Returns:
            object: The requested model instance, possibly enhanced with an imputation strategy if specified.

        Raises:
            KeyError: If the model with the specified name is not found in the library.
            ValueError: If an invalid imputation method is specified.

        Available imputation_method selections are:
            - "none": No imputation is applied.
            - "averaging": Imputes missing values by calculating the average of data in the input set.
            - "dunedin": Uses hybrid_impute with the dunedinPACE gold standard values.
            - "sesame_450k": Uses hybrid_impute with the sesame 450k gold standard values.
        If no imputation method is selected then the default imputation method specified within the model's definition is used.
        """
        if name not in self.model_definitions:
            raise KeyError(f"Model not found: {name}")

        model_def = self.model_definitions[name]
        model_def["name"] = name
        model_type = model_def["model"]["type"]
        model_instance = self.model_builders[model_type](model_def)

        global_default = "sesame_450k"
        default_imputation_method = model_def["model"].get(
            "default_imputation", global_default
        )
        imputation_method = imputation_method or default_imputation_method

        if imputation_method == "none":
            return model_instance
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
        elif imputation_method == "sesame_450k":
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
        for name, model_def in self.model_definitions.items():
            if (species is None or model_def.get("species") == species) and (
                tissue is None or model_def.get("tissue") == tissue
            ):
                matches[name] = model_def
        return matches
