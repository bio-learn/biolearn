"""
Test for HurdleAPIModel
"""

import pytest
import pandas as pd
import numpy as np
from unittest.mock import Mock, patch, MagicMock
import os

from biolearn.model import HurdleAPIModel
from biolearn.model_gallery import ModelGallery


class TestHurdleAPIModel:
    """Test the HurdleAPIModel class."""

    # Test CpG sites for mock data
    TEST_CPG_SITES = [
        "cg00000029",
        "cg00000321",
        "cg00000714",
        "cg00001793",
        "cg00002028",
        "cg00002426",
        "cg00002719",
        "cg00003091",
        "cg00003287",
        "cg00003994",
    ] + [f"cg{i:08d}" for i in range(10, 100)]

    def test_init_without_api_key(self):
        """Test that initialization fails without API key."""
        with pytest.raises(ValueError, match="API key required"):
            HurdleAPIModel()

    def test_init_with_api_key(self):
        """Test initialization with API key."""
        model = HurdleAPIModel(api_key="test_key")
        assert model.api_key == "test_key"
        assert "sandbox" in model.api_endpoint

    @patch.dict(os.environ, {"HURDLE_API_KEY": "env_test_key"})
    def test_init_with_env_key(self):
        """Test initialization with environment variable."""
        model = HurdleAPIModel()
        assert model.api_key == "env_test_key"

    def test_production_endpoint(self):
        """Test production endpoint selection."""
        model = HurdleAPIModel(api_key="test_key", use_production=True)
        assert "sandbox" not in model.api_endpoint
        assert "api.hurdle.bio" in model.api_endpoint

    @patch("builtins.input", return_value="no")
    def test_consent_denied(self, mock_input):
        """Test that prediction fails when consent is denied."""
        model = HurdleAPIModel(api_key="test_key")
        data = pd.DataFrame(np.random.rand(100, 5))

        with pytest.raises(ValueError, match="User consent required"):
            model.predict(data)

    def test_missing_cpgs_error(self):
        """Test that missing CpG sites raise an informative error."""
        model = HurdleAPIModel(api_key="test_key")

        # Create data missing some required CpGs (only use half of test CpGs)
        incomplete_cpgs = self.TEST_CPG_SITES[:50]
        data = pd.DataFrame(
            np.random.rand(50, 2),
            index=incomplete_cpgs,
            columns=["sample_0", "sample_1"],
        )

        # Mock loading required CpGs
        model.required_cpgs = self.TEST_CPG_SITES

        with pytest.raises(
            ValueError,
            match=r"Missing \d+/\d+ required CpG sites.*Please impute missing values",
        ):
            model.predict(data, require_consent=False)

    @patch("builtins.input", return_value="yes")
    @patch("requests.post")
    @patch("requests.put")
    def test_predict_success(self, mock_put, mock_post, mock_input):
        """Test successful prediction workflow."""
        # Setup model
        model = HurdleAPIModel(api_key="test_key")
        model.required_cpgs = self.TEST_CPG_SITES

        # Mock API responses
        mock_post.return_value.status_code = 200
        mock_post.return_value.json.side_effect = [
            {
                "uploadUrl": "https://example.com/upload",
                "requestId": "test-request-123",
            },
            {
                "data": [
                    {"barcode": "sample_0", "prediction": "35.5"},
                    {"barcode": "sample_1", "prediction": "42.3"},
                ]
            },
        ]

        mock_put.return_value.status_code = 200

        # Create test data with proper CpG indices
        data = pd.DataFrame(
            np.random.rand(100, 2),
            index=self.TEST_CPG_SITES,
            columns=["sample_0", "sample_1"],
        )

        # Make predictions
        predictions = model.predict(data)

        assert isinstance(predictions, pd.Series)
        assert len(predictions) == 2
        assert predictions["sample_0"] == 35.5
        assert predictions["sample_1"] == 42.3

    @patch("builtins.input", return_value="yes")
    @patch("requests.post")
    def test_api_error_handling(self, mock_post, mock_input):
        """Test API error handling."""
        model = HurdleAPIModel(api_key="test_key")
        model.required_cpgs = self.TEST_CPG_SITES

        # Mock failed API response
        mock_post.return_value.status_code = 401
        mock_post.return_value.text = "Unauthorized"

        # Create test data with proper CpG indices
        data = pd.DataFrame(
            np.random.rand(100, 2),
            index=self.TEST_CPG_SITES,
            columns=["sample_0", "sample_1"],
        )

        with pytest.raises(Exception, match="Failed to get upload URL: 401"):
            model.predict(data)

    def test_model_gallery_integration(self):
        """Test that HurdleAPIModel is available in ModelGallery."""
        gallery = ModelGallery()

        # Check that model is in definitions
        assert "HurdleInflammAge" in gallery.model_definitions

        # Test loading with API key
        with patch.dict(os.environ, {"HURDLE_API_KEY": "test_key"}):
            model = gallery.get("HurdleInflammAge")
            assert isinstance(model, HurdleAPIModel)

    @patch("builtins.input", return_value="yes")
    def test_consent_only_asked_once(self, mock_input):
        """Test that consent is only requested once."""
        model = HurdleAPIModel(api_key="test_key")
        model.required_cpgs = self.TEST_CPG_SITES

        # Mock successful API calls
        with (
            patch("requests.post") as mock_post,
            patch("requests.put") as mock_put,
        ):
            mock_post.return_value.status_code = 200
            mock_post.return_value.json.side_effect = [
                {"uploadUrl": "url", "requestId": "id"},
                {"data": [{"barcode": "sample", "prediction": "40"}]},
            ] * 2  # For two calls
            mock_put.return_value.status_code = 200

            # Create test data with proper CpG indices
            data = pd.DataFrame(
                np.random.rand(100, 1),
                index=self.TEST_CPG_SITES,
                columns=["sample"],
            )

            # First prediction
            model.predict(data)

            # Second prediction should not ask for consent again
            model.predict(data)

            # Input should only be called once
            assert mock_input.call_count == 1
