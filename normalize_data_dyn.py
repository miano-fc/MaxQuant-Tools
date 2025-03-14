# Import Modules
import numpy as np
import pandas as pd

# Read raw data and convert from .txt to .csv
def convert_txt(experimental_data_txt):
    print("Converting peptides.txt into a .csv")
    experimental_data = pd.read_csv(experimental_data_txt, delimiter="\t")
    experimental_data.to_csv("experimental_data/peptides.csv", index=False)
    return experimental_data

# User inputs dynamic number of experiment/control pairs
def get_experiments():
    experiment_list = []
    while True:
        experiment = input("Enter experiment name (or type 'done' to finish): ")
        if experiment.lower() == 'done':
            break
        control = input(f"Enter control name for {experiment}: ")
        experiment_list.append((experiment, control))
    print("\n")
    return experiment_list

# Remove NaN values from raw data and replace them with zeroes
def remove_nan(experiment_list, experimental_data):
    print('Replacing NaN values')
    for experiment, control in experiment_list:
        experimental_data[f'Experiment {experiment}'] = experimental_data[f'Experiment {experiment}'].fillna(0)
        experimental_data[f'Experiment {control}'] = experimental_data[f'Experiment {control}'].fillna(0)
    return experimental_data

# Subtract experimental 'Intensity' and 'Count' columns from control columns
def normalize_intensity(experiment_list, experimental_data):
    print('Normalizing intensity and count')
    experimental_data = remove_nan(experiment_list, experimental_data)
    
    for experiment, control in experiment_list:
        experimental_data[f'Count {experiment}'] = experimental_data[f'Experiment {experiment}'] - experimental_data[f'Experiment {control}']
        experimental_data[f'Intensity Experiment {experiment}'] = experimental_data[f'Intensity {experiment}'] - experimental_data[f'Intensity {control}']
    
    return experimental_data

# User selects significance and MS/MS count thresholds
def select_significance():
    return float(input("Enter PEP upper threshold: "))

def select_msms_count():
    return float(input("Enter MS/MS count lower threshold: "))

# Apply filtering criteria
def filter_data(experiment_list, experimental_data):
    print('Filtering data')
    experimental_data = normalize_intensity(experiment_list, experimental_data)
    print('\n')
    significance_threshold = select_significance()
    msms_count_threshold = select_msms_count()
    print('\n')
    
    experimental_data = experimental_data[
        (experimental_data['MS/MS Count'] >= msms_count_threshold) &
        (experimental_data['PEP'] < significance_threshold) &
        (experimental_data['Potential contaminant'].isna())
    ]
    
    # Remove rows only if ALL experiments have Intensity ≤ 0 and Count = 0
    mask_intensity = ~(experimental_data[[f'Intensity {exp}' for exp, _ in experiment_list]] <= 0).all(axis=1)
    mask_count = ~(experimental_data[[f'Count {exp}' for exp, _ in experiment_list]] == 0).all(axis=1)
    experimental_data = experimental_data.loc[mask_intensity & mask_count]
    
    return experimental_data

# Combine rows containing the same protein
def combine_rows(experiment_list, experimental_data):
    print('Aggregating data')
    experimental_data = filter_data(experiment_list, experimental_data)
    grouped_data = experimental_data.groupby('Protein names').agg(lambda x: x.sum() if np.issubdtype(x.dtype, np.number) else x.iloc[0])
    grouped_data.reset_index(inplace=True)
    return grouped_data

# Remove extra columns dynamically
def remove_extra_columns(experiment_list, experimental_data):
    print('Removing extra columns')
    experimental_data = combine_rows(experiment_list, experimental_data)
    
    # List of general columns to drop
    columns_to_drop = [
        'Sequence', 'N-term cleavage window', 'C-term cleavage window',
        'Amino acid before', 'First amino acid', 'Second amino acid',
        'Second last amino acid', 'Last amino acid', 'Amino acid after',
        'Length', 'Missed cleavages', 'Mass', 'Leading razor protein',
        'Start position', 'End position', 'Charges', 'PEP', 'Score',
        'Reverse', 'Potential contaminant', 'id', 'Protein group IDs',
        'Mod. peptide IDs', 'Evidence IDs', 'MS/MS IDs', 'Best MS/MS',
        'Taxonomy IDs', 'Taxonomy names', 'Mass deficit', 'Unique (Groups)',
        'Unique (Proteins)', 'Deamidation (N) site IDs', 'A Count', 'R Count',
        'N Count', 'D Count', 'C Count', 'Q Count', 'E Count', 'G Count',
        'H Count', 'I Count', 'L Count', 'K Count', 'M Count', 'F Count',
        'P Count', 'S Count', 'T Count', 'W Count', 'Y Count', 'V Count',
        'U Count', 'O Count',

    ]
    
    # Dynamically remove experiment-specific columns
    for experiment, control in experiment_list:
        columns_to_drop.extend([
            f'Intensity {experiment}', f'Intensity {control}',
            f'Experiment {experiment}', f'Experiment {control}'
        ])
    
    experimental_data = experimental_data.loc[:, ~experimental_data.columns.isin(columns_to_drop)]
    experimental_data.to_csv("results/aggregated_data.csv", index=False)
    return experimental_data

# Define a main function
def main(experimental_data_txt):
    print("Normalize MaxQuant peptides.txt by Intensity and aggregate by Protein names\n")
    experiment_list = get_experiments()
    experimental_data = convert_txt(experimental_data_txt)
    remove_extra_columns(experiment_list, experimental_data)
    print("Done.")
    return remove_extra_columns
