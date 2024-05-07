# This script is used together with the "py_env" environment
import numpy as np
import umap
from scipy.spatial.distance import pdist, squareform
import pandas as pd
import glob
import os

# Directory path where CSV files with abundance data are stored
input_directory_path = '/PATH/umap_input_data/some_directory/'

# Directory path where you want to save the output CSV files
output_directory_path = '/PATH/umap_embedding_data/some_directory/' # Change this when altering parameters

# Ensure the output directory exists
if not os.path.exists(output_directory_path):
    os.makedirs(output_directory_path)
    print("Created output directory")

# Automatically list all CSV files in the input directory
file_paths = glob.glob(input_directory_path + '*.csv')

# Loop over each file path
for file_path in file_paths:
    # Read the data
    data = pd.read_csv(file_path)
    
    # Compute Bray-Curtis distance matrix
    bray_curtis_matrix = squareform(pdist(data, metric='braycurtis'))
    
    # Defining UMAP parameters
    n_components = 10 # Default = 2
    min_dist = 0.1   # Default = 0.1
    n_neighbors = 250 # Default = 15
    random_state = 123
    
    # Apply UMAP with the precomputed distance matrix
    reducer = umap.UMAP(metric='precomputed', n_components=n_components, min_dist=min_dist, n_neighbors=n_neighbors, random_state=random_state)
    embedding = reducer.fit_transform(bray_curtis_matrix)
    
    # Create list of column names
    column_names = ['UMAP_{}'.format(i) for i in range(1, n_components + 1)]
    
    # Convert embedding to DataFrame
    embedding_df = pd.DataFrame(embedding, columns=column_names)
    
    # Extract the base name of the file without the directory path
    base_file_name = os.path.basename(file_path)
    
    # Full path for the output file
    full_output_path = os.path.join(output_directory_path, base_file_name)
    
    # Save to CSV
    embedding_df.to_csv(full_output_path, index=False)
