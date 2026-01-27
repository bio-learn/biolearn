# Using Hurdle's InflammAge Model in Biolearn

## Non-Commercial Use Only

**Important Notice:**
This inflammaging calculation is provided for personal and educational purposes only.

Any commercial use, distribution, or resale of this information is strictly prohibited without prior written agreement from Chronomics Limited. See https://hurdle.bio/ for further commercial information.

---

## Quick Start

1. **Get API Credentials**
   - Register at: https://dashboard.hurdle.bio/register
   - Navigate to Developer > API keys
   - Save your API key securely

2. **Set Up Your API Key**
   ```bash
   export HURDLE_API_KEY="your_api_key_here"
   ```

3. **Use the Model**
   ```python
   from biolearn.model_gallery import ModelGallery
   
   gallery = ModelGallery()
   model = gallery.get("HurdleInflammAge")
   
   # Make predictions
   predictions = model.predict(methylation_data)
   ```

## Privacy Notice

⚠️ **Data Privacy Warning**: Using this model will send your methylation data to Hurdle Bio's servers. You will be asked for explicit consent before any data is transmitted. Ensure you have permission to share genomic data with third parties.

## Detailed Setup

### Getting Your API Key

1. Visit https://dashboard.hurdle.bio/register
2. Create an account
3. Navigate to Developer > API keys
4. Generate a new API key
5. Store it securely (never commit to version control)

### Setting Your API Key

**Option 1: Environment Variable (Recommended)**
```bash
export HURDLE_API_KEY="your_api_key_here"
```

**Option 2: Direct Input**
```python
from biolearn.model import HurdleAPIModel

model = HurdleAPIModel(api_key="your_api_key_here")
```

## Usage Examples

### Basic Usage
```python
from biolearn.model_gallery import ModelGallery
from biolearn.data_library import DataLibrary

# Load data
data = DataLibrary().get("YourDataset").load()

# Get model
gallery = ModelGallery()
model = gallery.get("HurdleInflammAge")

# Make predictions (you'll be asked for consent)
predictions = model.predict(data.dnam)
```

### With Metadata
```python
# If your data includes age and sex information,
# the model will use this for more accurate predictions
predictions = model.predict(data)  # Uses both dnam and metadata
```

## Data Requirements

- **Input**: DNA methylation beta values (0-1 range)
- **Format**: Pandas DataFrame with CpG sites as rows, samples as columns
- **Tissue**: Blood samples recommended
- **CpG Sites**: Requires all 62,000+ specific CpG sites (no missing values allowed)
- **Optional**: Age and sex metadata improve accuracy

### Missing Value Handling

The model requires approximately 62,000 specific CpG sites. If your data is missing any of these sites, you'll get an error listing the missing CpGs. You must impute missing values before calling the model. Use biolearn's built-in imputation:

```python
# Option 1: Use ModelGallery's imputation
model = gallery.get("HurdleInflammAge", imputation_method="averaging")

# Option 2: Impute your data beforehand
from biolearn.imputation import hybrid_impute
imputed_data = hybrid_impute(your_data, reference_values, required_cpgs)
```

## Troubleshooting

### Common Errors

**"API key required"**
- Set your API key using one of the methods above
- Verify the key is correct

**"Failed to get upload URL: 401"**
- Your API key may be invalid or expired
- Verify you're using a production API key from https://dashboard.hurdle.bio

**"Missing X/Y required CpG sites"**
- Your data is missing some required CpG sites
- The error will show examples of missing sites
- Use one of biolearn's imputation methods before calling predict()

**"No required CpG sites found in data"**
- Your methylation data doesn't contain any of the required CpG sites  
- Ensure your data uses Illumina array probe IDs (e.g., "cg00000029")

**"User consent required to send data to external API"**
- You must type "yes" when prompted to consent to sending data
- This ensures you're aware that data will be transmitted to external servers

**Network Errors**
- Check your internet connection
- Verify firewall settings allow HTTPS connections to hurdle.bio

## Best Practices

1. **Test First**: Verify your API credentials work with a small subset of samples
2. **Batch Processing**: For large datasets, consider processing in smaller batches
3. **Error Handling**: Always wrap API calls in try-except blocks
4. **Data Security**: Ensure you have consent to share genomic data before using this model
5. **Save Results**: Cache predictions locally to avoid redundant API calls

## Support

- **Biolearn Issues**: Open an issue on the [biolearn GitHub repository](https://github.com/bio-learn/biolearn)
- **API Access**: Contact Hurdle Bio support via https://hurdle.bio
- **Commercial Use**: Contact Chronomics Limited for licensing inquiries

## Example Script

See `examples/hurdle_api_example.py` for a complete working example.

## Technical Details

- **Model Type**: Third-party API integration
- **Output**: InflammAge score (inflammation-based biological age)
- **API Endpoint**: https://api.hurdle.bio/predict/v1/
- **Timeout**: 30 seconds default
- **Consent**: Required once per session