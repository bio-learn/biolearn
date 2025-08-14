# Using Hurdle's Inflammage Model in Biolearn

## Quick Start

1. **Get API Credentials**
   - Register at: https://dashboard.sandbox.hurdle.bio/register/partner
   - Save your API key securely

2. **Set Up Your API Key**
   ```bash
   export HURDLE_API_KEY="your_api_key_here"
   ```

3. **Use the Model**
   ```python
   from biolearn.model_gallery import ModelGallery
   
   gallery = ModelGallery()
   model = gallery.get("HurdleInflammage")
   
   # Make predictions
   predictions = model.predict(methylation_data)
   ```

## Important Privacy Notice

⚠️ **Data Privacy Warning**: Using this model will send your methylation data to Hurdle Bio's servers. You will be asked for explicit consent before any data is transmitted. Ensure you have permission to share genomic data with third parties.

## Detailed Setup

### Getting Your API Key

1. Visit https://dashboard.sandbox.hurdle.bio/register/partner
2. Create an account
3. Navigate to API Keys section
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
model = gallery.get("HurdleInflammage")

# Make predictions (you'll be asked for consent)
predictions = model.predict(data.dnam)
```

### With Metadata
```python
# If your data includes age and sex information
# The model will use this for more accurate predictions
predictions = model.predict(data)  # Uses both dnam and metadata
```

## Troubleshooting

### Common Errors

**"API key required"**
- Set your API key using one of the methods above
- Verify the key is correct

**"Failed to get upload URL: 401"**
- Your API key may be invalid or expired
- Check you're using the correct environment (sandbox vs production)

**"No required CpG sites found in data"**
- Your methylation data doesn't contain the CpG sites Hurdle requires
- Contact Hurdle Bio for their required CpG sites list

**Network Errors**
- Check your internet connection
- Verify firewall settings allow HTTPS connections to hurdle.bio

## Data Requirements

- **Input**: DNA methylation beta values (0-1 range)
- **Format**: Pandas DataFrame with CpG sites as rows, samples as columns
- **Tissue**: Blood samples recommended
- **Optional**: Age and sex metadata improve accuracy

## Sandbox vs Production

By default, the model uses Hurdle's sandbox environment for testing:

```python
# Sandbox (default)
model = gallery.get("HurdleInflammage")

# Production (requires production API key)
model = HurdleAPIModel(api_key="prod_key", use_production=True)
```

## Best Practices

1. **Test First**: Use sandbox environment before production
2. **Batch Processing**: For large datasets, process in smaller batches
3. **Error Handling**: Always wrap API calls in try-except blocks
4. **Data Security**: Ensure you have consent to share genomic data

## Support

- **Biolearn Issues**: Open an issue on the biolearn GitHub repository
- **API Access**: Contact Hurdle Bio support
- **Commercial Use**: Contact Hurdle Bio for licensing

## Example Script

See `examples/hurdle_api_example.py` for a complete working example.