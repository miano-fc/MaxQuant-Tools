import numpy as np
import pandas as pd
import normalize_data_dyn as norm
import normalize_data_no_controls_dyn as norm_exponly
import sort_data_dyn as sort
import common_data as common

# Ask whether the user has experiments only or experiments with corresponding controls
def get_data_type():
    data_type = int(input("Enter your data type (1/2): \n1. Experiments only \n2. Experiments with corresponding controls\n "))
    return data_type

# Process experiemntal data based on get_data_type_response and ouput aggregated_data.csv to the results folder
def preprocess_data(experimental_data):
    data_type = get_data_type()

    if data_type == 1:
        norm_exponly.main(experimental_data)
        
    if data_type == 2:
        # If a user's experiments have corresponding controls, normalize experimental data to control data
        norm.main(experimental_data)
    
    return experimental_data

raw_data = "experimental_data/peptides.txt" # Path to your peptides.txt file
preprocess_data(raw_data)

# Sort aggregated_data.csv by experiment/differential expression
# aggregated_data = pd.read_csv("results/aggregated_data.csv")
# sort.main(aggregated_data)
