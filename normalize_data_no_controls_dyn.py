import numpy as np
import pandas as pd
import normalize_data_dyn as norm
import common_data as common

# Script for cleaning MaxQuant data containing experiments w/o corresponding controls
# 1. Get a list of experiments from the user
def get_experiments():
    experiment_list = []
    
    while True:
        experiment = input("Enter experiment name or enter 'done' to finish: ")
        if experiment.lower() == 'done':
            break

        experiment_list.append(experiment)
    print("\n")

    return experiment_list

# 2. Remove NaN values from peptides.txt and replace them wih zeroes
def remove_nan_values(experiment_list, experimental_data):
    print("Removing NaN values")

    for experiment in experiment_list:
        experimental_data[f'Experiment {experiment}'] = experimental_data[f'Experiment {experiment}'].fillna(0)
    
    return experimental_data

# 3. User selects MS/MS count and PEP significance thresholds and rows where 'Potential contminant' = '+' are removed
def filter_data(experiment_list, experimental_data):
    print("Filtering experimental data")
    experimental_data = remove_nan_values(experiment_list, experimental_data)

    print("\n")    
    significance_threshold = norm.select_significance()
    msms_count_threshold = norm.select_msms_count()
    print("\n")

    experimental_data = experimental_data[
        (experimental_data['MS/MS Count'] >= msms_count_threshold) &
        (experimental_data['PEP'] <= significance_threshold) &
        (experimental_data['Potential contaminant'].isna())
    ]

    # Remove rows only if ALL experiments have Intensity ≤ 0 and Count = 0
    mask_intensity = ~(experimental_data[[f'Intensity {exp}' for exp in experiment_list]] <= 0).all(axis=1)
    mask_count = ~(experimental_data[[f'Experiment {exp}' for exp in experiment_list]] == 0).all(axis=1)
    experimental_data = experimental_data.loc[mask_intensity & mask_count]

    return experimental_data

# 4. Combine rows containing the same 'Protein name' value
def group_rows(experiment_list, experimental_data):
    print("Aggregating rows by protein name")
    experimental_data = filter_data(experiment_list, experimental_data)
    grouped_data = experimental_data.groupby('Protein names').agg(lambda x: x.sum() if np.issubdtype(x.dtype, np.number) else x.iloc[0])
    grouped_data.reset_index(inplace=True)

    return grouped_data

# 5. Drop unnecessary columns
def drop_unnecessary_columns(experiment_list, experimental_data):
    print("Dropping unnecessary columns")
    experimental_data = group_rows(experiment_list, experimental_data)
    columns_to_drop = common.columns_to_drop

    experimental_data = experimental_data.loc[:, ~experimental_data.columns.isin(columns_to_drop)]
    experimental_data.to_csv("results/aggregated_data.csv", index=False)
    return experimental_data

# 6. Define a main function
def main(experimental_data_txt):
    print("Filter MaxQuant data and aggregate rows by protein name")

    experiment_list = get_experiments()
    experimental_data = norm.convert_txt(experimental_data_txt)

    drop_unnecessary_columns(experiment_list, experimental_data)
    print("Done.")

    return experimental_data


