# Import Modules
import numpy as np
import pandas as pd

# Read raw data
def get_experiments():
    experiment_1 = input("Enter experiment 1: ")
    control_1 = input("Enter control 1: ")
    experiment_2 = input("Enter experiment 2: ")
    control_2 = input("Enter control 2: ")

    experiment_list = [experiment_1, control_1, experiment_2, control_2]

    return experiment_list

# Remove NaN values from raw data and replace them with zeroes
def remove_nan(experiment_list, experimental_data):
    print('Replacing NaN values')
    for experiment in experiment_list:
        experimental_data['Experiment ' + experiment] = experimental_data['Experiment ' + experiment].fillna(0)

    return experimental_data

# Subtract experimental 'Intensity' and 'Count' columns from control columns
def normalize_intensity(experiment_list, experimental_data):
    print('Normalizing intensity and count')
    experimental_data = remove_nan(experiment_list, experimental_data)
    
    experimental_data.loc[:, 'Count ' + experiment_list[0]] = experimental_data.loc[:, 'Experiment ' + experiment_list[0]] - experimental_data.loc[:, 'Experiment ' + experiment_list[1]]
    experimental_data.loc[:, 'Count ' + experiment_list[2]] = experimental_data.loc[:, 'Experiment ' + experiment_list[2]] - experimental_data.loc[:, 'Experiment ' + experiment_list[3]]

    experimental_data.loc[:, 'Intensity Experiment ' + experiment_list[0]] = experimental_data.loc[:, 'Intensity ' + experiment_list[0]] - experimental_data.loc[:, 'Intensity ' + experiment_list[1]]
    experimental_data.loc[:, 'Intensity Experiment ' + experiment_list[2]] = experimental_data.loc[:, 'Intensity ' + experiment_list[2]] - experimental_data.loc[:, 'Intensity ' + experiment_list[3]]

    return experimental_data

# Remove rows that meet the following criteria
# 1. Contaminant-containing
# 2. PEP > user-selected threshold (insignificant values)
# 3. MS/MS count >= user-selected threshold
# 4. Experiment 1 and 2 intensity are both 0 or negative
# 5. Experiment 1 and 2 count are both 0
def select_significance():
    significance = float(input("Enter PEP upper threshold: "))
    return significance

def select_msms_count():
    msms_count = float(input("Enter MS/MS count lower threshold: "))
    return msms_count

def filter_data(experiment_list, experimental_data):
    print('Filtering data')
    experimental_data = normalize_intensity(experiment_list, experimental_data)

    significance_threshold = select_significance()
    msms_count_threshold = select_msms_count()

    experimental_data = experimental_data[
        (experimental_data['MS/MS Count'] >= msms_count_threshold) &
        (experimental_data['PEP'] < significance_threshold) &
        (experimental_data['Potential contaminant'].isna())
        ]

    mask_intensity = ~((experimental_data['Intensity ' + experiment_list[0]] <= 0) & (experimental_data['Intensity ' + experiment_list[2]] <= 0))
    experimental_data = experimental_data.loc[mask_intensity]

    mask_count = ~((experimental_data['Count ' + experiment_list[0]] == 0) & (experimental_data['Count ' + experiment_list[2]] == 0))
    experimental_data = experimental_data[mask_count]  # Keep rows where at least one count is nonzero

    return experimental_data

# Combine rows containing the same protein
def combine_rows(experiment_list, experimental_data):
    print('Aggregating data')
    experimental_data = filter_data(experiment_list, experimental_data)

    grouped_data = experimental_data.groupby('Protein names')
    aggregated_data = grouped_data.agg(lambda x: x.sum() if np.issubdtype(x.dtype, np.number) else x.iloc[0])
    aggregated_data.reset_index(inplace=True)

    return aggregated_data

# Remove extra columns
def remove_extra_columns(experiment_list, experimental_data):
    print('Removing extra columns')
    experimental_data = combine_rows(experiment_list, experimental_data)

    columns_to_drop = ['Sequence',
            'N-term cleavage window',
            'C-term cleavage window',
            'Amino acid before',
            'First amino acid',
            'Second amino acid',
            'Second last amino acid',
            'Last amino acid',
            'Amino acid after',
            'A Count',
            'R Count',
            'N Count',
            'D Count',
            'C Count',
            'Q Count',
            'E Count',
            'G Count',
            'H Count',
            'I Count',
            'L Count',
            'K Count',
            'M Count',
            'F Count',
            'P Count',
            'S Count',
            'T Count',
            'W Count',
            'Y Count',
            'V Count',
            'U Count',
            'O Count',
            'Length',
            'Missed cleavages',
            'Mass',
            'Leading razor protein',
            'Start position',
            'End position',
            'Charges',
            'PEP',
            'Score',
            'Reverse',
            'Potential contaminant',
            'id',
            'Protein group IDs',
            'Mod. peptide IDs',
            'Evidence IDs',
            'MS/MS IDs',
            'Best MS/MS',
            'Oxidation (M) site IDs',
            'Taxonomy IDs',
            'Taxonomy names',
            'Mass deficit',
            'Unique (Groups)',
            'Unique (Proteins)',
            'Deamidation (N) site IDs',
            'Intensity ' + experiment_list[0],
            'Intensity ' + experiment_list[1],
            'Intensity ' + experiment_list[2],
            'Intensity ' + experiment_list[3],
            'Experiment ' + experiment_list[0],
            'Experiment ' + experiment_list[1],
            'Experiment ' + experiment_list[2],
            'Experiment ' + experiment_list[3]
            ]

    experimental_data = experimental_data.loc[:, ~experimental_data.columns.isin(columns_to_drop)]
    experimental_data.to_csv("results/aggregated_data.csv", index = False)

    return experimental_data

# Define a main function
def main(experimental_data):
    experiment_list = get_experiments()

    remove_extra_columns(experiment_list, experimental_data)

    print("Done.")

    return remove_extra_columns

