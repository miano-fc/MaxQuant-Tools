import numpy as np
import pandas as pd
import normalize_raw_data_options as norm
import sort_experimental_data as sort

# Read raw data, normalize by experiment 1/2 count and intensity, then output aggregated_data.csv to the results folder
raw_data = "experimental_data/peptides.txt" # Path to your peptides.txt file
norm.main(raw_data)

# Sort aggregated_data.csv by experiment/differential expression
aggregated_data = pd.read_csv("results/aggregated_data.csv")
sort.main(aggregated_data)
