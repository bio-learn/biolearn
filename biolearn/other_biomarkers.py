from biolearn.util import get_data_file
import pandas as pd
import numpy as np

def estimate_sex(dnam_data):
    # Load the reference data
    reference = pd.read_csv(get_data_file("estimateSex.csv"), index_col=0)
    
    # Ensure the input data is a pandas DataFrame
    if not isinstance(dnam_data, pd.DataFrame):
        raise ValueError("Input DNAm data must be a pandas DataFrame")

    # Filter the reference data and input data to have matching probe IDs
    common_probes = dnam_data.index.intersection(reference.index)
    reference = reference.loc[common_probes]
    dnam_data = dnam_data.loc[common_probes]
    
    # Calculate the mean and standard deviation for autosomes
    autosomes = reference.loc[~reference['CHR'].isin(['X', 'Y'])].index
    auto_data = dnam_data.loc[autosomes]
    d_mean = auto_data.mean(axis=0, skipna=True)
    d_std = auto_data.std(axis=0, skipna=True)
    
    # Normalize the β-values using Z-score normalization
    z_data = (dnam_data.subtract(d_mean, axis=1)).div(d_std, axis=1)
    
    # Perform sex prediction using PCA coefficients from the reference data
    pred_xy = {}
    for chr in ['X', 'Y']:
        chr_ref = reference.loc[(reference['pca'] == chr)]
        chr_data = z_data.loc[chr_ref.index]
        chr_data.fillna(0, inplace=True)
        
        pred = (chr_data.T - chr_ref['mean'].values).dot(chr_ref['coeff'].values)
        pred_xy[chr] = pred
    
    # Create a DataFrame with the predictions
    pred_df = pd.DataFrame(pred_xy)
    
    # Determine the predicted sex based on the X and Y predictions
    pred_df['predicted_sex'] = 'Female'
    pred_df.loc[(pred_df['X'] < 0) & (pred_df['Y'] > 0), 'predicted_sex'] = 'Male'
    pred_df.loc[(pred_df['X'] > 0) & (pred_df['Y'] > 0), 'predicted_sex'] = '47,XXY'
    pred_df.loc[(pred_df['X'] < 0) & (pred_df['Y'] < 0), 'predicted_sex'] = '45,XO'
    
    return pred_df