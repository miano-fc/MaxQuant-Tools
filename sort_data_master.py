#!/usr/bin/env python3
"""
MaxQuant Proteomics Data Analysis Master Program

"""

import os
import numpy as np
import pandas as pd
import argparse

# Create directories if they don't exist
os.makedirs("experimental_data", exist_ok=True)
os.makedirs("results", exist_ok=True)

# ===== DATA CONVERSION AND PREPARATION FUNCTIONS =====

def convert_txt(experimental_data_txt):
    """Convert peptides.txt to CSV format"""
    print("Converting peptides.txt into a .csv")
    experimental_data = pd.read_csv(experimental_data_txt, delimiter="\t")
    experimental_data.to_csv("experimental_data/peptides.csv", index=False)
    return experimental_data

# ===== USER INPUT FUNCTIONS =====

def get_analysis_type():
    """Ask the user which type of analysis to perform"""
    print("\nSelect analysis type:")
    print("1: Experiments with corresponding controls")
    print("2: Experiments without controls")
    
    while True:
        choice = input("Enter your choice (1 or 2): ")
        if choice in ['1', '2']:
            return int(choice)
        print("Invalid choice. Please enter 1 or 2.")

def get_experiments_with_controls():
    """Get experiment/control pairs from user"""
    experiment_list = []
    print("\nEnter experiment/control pairs (type 'done' to finish):")
    while True:
        experiment = input("Enter experiment name (or type 'done' to finish): ")
        if experiment.lower() == 'done':
            break
        control = input(f"Enter control name for {experiment}: ")
        experiment_list.append((experiment, control))
    print("\n")
    return experiment_list

def get_experiments_without_controls():
    """Get experiment names from user without controls"""
    experiment_list = []
    print("\nEnter experiment names (type 'done' to finish):")
    while True:
        experiment = input("Enter experiment name (or type 'done' to finish): ")
        if experiment.lower() == 'done':
            break
        experiment_list.append(experiment)
    print("\n")
    return experiment_list

def select_significance():
    """Get PEP threshold from user"""
    return float(input("Enter PEP upper threshold (e.g. 0.05): "))

def select_msms_count():
    """Get MS/MS count threshold from user"""
    return float(input("Enter MS/MS count lower threshold (e.g. 2): "))

def select_log2fc(experiment_list):
    """Let user select how Log2 fold change should be calculated"""
    print("\nSelect how Log2 fold change should be calculated:")
    print("Options:")
    for i in range(len(experiment_list)):
        for j in range(i + 1, len(experiment_list)):
            print(f"  {i + 1}/{j + 1} ( {experiment_list[i]} / {experiment_list[j]} )")
            print(f"  {j + 1}/{i + 1} ( {experiment_list[j]} / {experiment_list[i]} )")
    
    selection = input("Enter your choice (e.g., 1/2 or 2/1): ")
    return selection

# ===== DATA PROCESSING FUNCTIONS =====

def remove_nan_with_controls(experiment_list, experimental_data):
    """Replace NaN values with zeros for experiments with controls"""
    print('Replacing NaN values')
    
    for experiment, control in experiment_list:
        experimental_data[f'Experiment {experiment}'] = experimental_data[f'Experiment {experiment}'].fillna(0)
        experimental_data[f'Experiment {control}'] = experimental_data[f'Experiment {control}'].fillna(0)
    return experimental_data

def remove_nan_without_controls(experiment_list, experimental_data):
    """Replace NaN values with zeros for experiments without controls"""
    print("Removing NaN values")
    for experiment in experiment_list:
        experimental_data[f'Experiment {experiment}'] = experimental_data[f'Experiment {experiment}'].fillna(0)
    
    return experimental_data

def normalize_intensity(experiment_list, experimental_data):
    """Subtract control from experiment intensities"""
    print('Normalizing intensity and count')
    experimental_data = remove_nan_with_controls(experiment_list, experimental_data)
    
    for experiment, control in experiment_list:
        experimental_data[f'Count {experiment}'] = experimental_data[f'Experiment {experiment}'] - experimental_data[f'Experiment {control}']
        experimental_data[f'Intensity Experiment {experiment}'] = experimental_data[f'Intensity {experiment}'] - experimental_data[f'Intensity {control}']
    
    return experimental_data

def filter_data_with_controls(experiment_list, experimental_data):
    """Filter data for experiments with controls"""
    print('Filtering data')
    experimental_data = normalize_intensity(experiment_list, experimental_data)
    print('\n')
    significance_threshold = select_significance()
    msms_count_threshold = select_msms_count()
    print('\n')
    
    experimental_data = experimental_data[
        (experimental_data['MS/MS Count'] >= msms_count_threshold) &
        (experimental_data['PEP'] <= significance_threshold) &
        (experimental_data['Potential contaminant'].isna())
    ]
    
    # Remove rows only if ALL experiments have Intensity ≤ 0 and Count = 0
    mask_intensity = ~(experimental_data[[f'Intensity {exp}' for exp, _ in experiment_list]] <= 0).all(axis=1)
    mask_count = ~(experimental_data[[f'Count {exp}' for exp, _ in experiment_list]] == 0).all(axis=1)
    experimental_data = experimental_data.loc[mask_intensity & mask_count]
    
    return experimental_data

def filter_data_without_controls(experiment_list, experimental_data):
    """Filter data for experiments without controls"""
    print("Filtering experimental data")
    experimental_data = remove_nan_without_controls(experiment_list, experimental_data)
    print("\n")    
    significance_threshold = select_significance()
    msms_count_threshold = select_msms_count()
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

def combine_rows_with_controls(experiment_list, experimental_data):
    """Aggregate rows by protein name for experiments with controls"""
    print('Aggregating rows by protein name')
    experimental_data = filter_data_with_controls(experiment_list, experimental_data)
    grouped_data = experimental_data.groupby('Protein names').agg(lambda x: x.sum() if np.issubdtype(x.dtype, np.number) else x.iloc[0])
    grouped_data.reset_index(inplace=True)

    return grouped_data

def combine_rows_without_controls(experiment_list, experimental_data):
    """Aggregate rows by protein name for experiments without controls"""
    print("Aggregating rows by protein name")
    experimental_data = filter_data_without_controls(experiment_list, experimental_data)
    grouped_data = experimental_data.groupby('Protein names').agg(lambda x: x.sum() if np.issubdtype(x.dtype, np.number) else x.iloc[0])
    grouped_data.reset_index(inplace=True)
    return grouped_data

def remove_extra_columns_with_controls(experiment_list, experimental_data):
    """Remove unnecessary columns for experiments with controls"""
    print('Removing extra columns')
    experimental_data = combine_rows_with_controls(experiment_list, experimental_data)
    
    # Define columns to drop
    columns_to_drop = [
        'Peptide IDs', 'Peptide is razor', 'Mod. peptide IDs', 'Evidence IDs', 
        'MS/MS IDs', 'Best MS/MS', 'Reverse', 'Potential contaminant', 'id', 
        'Peptide counts (all)', 'Peptide counts (razor unique)', 'Protein IDs',
        'Peptide counts (unique)', 'Sequence coverage [%]', 'Unique + razor sequence coverage [%]',
        'Sequence coverage [%]', 'Mol. weight [kDa]', 'Sequence length', 'Q-value',
        'Score', 'Fasta headers'
    ]
    
    # Dynamically remove experiment-specific columns
    for experiment, control in experiment_list:
        columns_to_drop.extend([
            f'Intensity {experiment}', f'Intensity {control}',
            f'Experiment {experiment}', f'Experiment {control}'
        ])
    
    # Filter columns to drop only those that actually exist in the dataframe
    columns_to_drop = [col for col in columns_to_drop if col in experimental_data.columns]
    
    experimental_data = experimental_data.drop(columns=columns_to_drop, errors='ignore')
    experimental_data.to_csv("results/aggregated_data.csv", index=False)
    return experimental_data

def remove_extra_columns_without_controls(experiment_list, experimental_data):
    """Remove unnecessary columns for experiments without controls"""
    print("Dropping unnecessary columns")
    experimental_data = combine_rows_without_controls(experiment_list, experimental_data)
    
    # Define columns to drop
    columns_to_drop = [
        'Peptide IDs', 'Peptide is razor', 'Mod. peptide IDs', 'Evidence IDs', 
        'MS/MS IDs', 'Best MS/MS', 'Reverse', 'Potential contaminant', 'id', 
        'Peptide counts (all)', 'Peptide counts (razor unique)', 'Protein IDs',
        'Peptide counts (unique)', 'Sequence coverage [%]', 'Unique + razor sequence coverage [%]',
        'Sequence coverage [%]', 'Mol. weight [kDa]', 'Sequence length', 'Q-value',
        'Score', 'Fasta headers'
    ]
    
    # Filter columns to drop only those that actually exist in the dataframe
    columns_to_drop = [col for col in columns_to_drop if col in experimental_data.columns]
    
    experimental_data = experimental_data.drop(columns=columns_to_drop, errors='ignore')
    experimental_data.to_csv("results/aggregated_data.csv", index=False)
    return experimental_data

# ===== ANALYSIS FUNCTIONS =====

def get_differentially_expressed(experiment_list, experimental_data):
    """Generate a list of differentially expressed proteins"""
    print("Generating a list of differentially expressed proteins.\n")
    
    # For with_controls case, experiment_list contains tuples (exp, control)
    # For without_controls case, experiment_list contains experiment names only
    # Handle both cases
    if isinstance(experiment_list[0], tuple):
        # Extract experiment names from experiment_list tuples
        exp_names = [exp for exp, _ in experiment_list]
        
        condition = (experimental_data[[f'Count {exp}' for exp in exp_names]] != 0).all(axis=1) & \
                    (experimental_data[[f'Intensity Experiment {exp}' for exp in exp_names]] > 0).all(axis=1)
    else:
        # For without_controls, we need to create intensity columns first
        for exp in experiment_list:
            experimental_data[f'Intensity Experiment {exp}'] = experimental_data[f'Intensity {exp}']
            experimental_data[f'Count {exp}'] = experimental_data[f'Experiment {exp}']
        
        condition = (experimental_data[[f'Count {exp}' for exp in experiment_list]] != 0).all(axis=1) & \
                    (experimental_data[[f'Intensity Experiment {exp}' for exp in experiment_list]] > 0).all(axis=1)
        
        exp_names = experiment_list
    
    differential_expression = experimental_data.loc[condition]
    
    if len(exp_names) < 2:
        print("Need at least two experiments to calculate fold change.")
        return differential_expression
    
    user_selection = select_log2fc(exp_names)
    exp1_idx, exp2_idx = map(int, user_selection.split('/'))
    exp1, exp2 = exp_names[exp1_idx - 1], exp_names[exp2_idx - 1]
    
    differential_expression[f'Log2FC {exp1}/{exp2}'] = np.log2(
        differential_expression[f'Intensity Experiment {exp1}'] /
        differential_expression[f'Intensity Experiment {exp2}']
    )
    
    differential_expression.to_csv("results/differential_expression.csv", index=False)
    return differential_expression

def get_experiment_exclusive(experiment_list, experimental_data, exp_index):
    """Generate a list of proteins expressed exclusively in one experiment"""
    # Handle both with_controls and without_controls cases
    if isinstance(experiment_list[0], tuple):
        exp = experiment_list[exp_index][0]  # Extract experiment name from the tuple
        exp_names = [exp for exp, _ in experiment_list]
    else:
        exp = experiment_list[exp_index]
        exp_names = experiment_list
    
    print(f"Generating a list of proteins expressed exclusively in {exp}.")
    
    # Create intensity columns if they don't exist (for without_controls case)
    if f'Count {exp}' not in experimental_data.columns:
        for e in exp_names:
            experimental_data[f'Intensity Experiment {e}'] = experimental_data[f'Intensity {e}']
            experimental_data[f'Count {e}'] = experimental_data[f'Experiment {e}']
    
    condition = (experimental_data[f'Count {exp}'] > 0) & \
                (experimental_data[f'Intensity Experiment {exp}'] > 0)
    
    for i, other_exp in enumerate(exp_names):
        if other_exp != exp:
            condition &= (experimental_data[f'Count {other_exp}'] <= 0) & \
                         (experimental_data[f'Intensity Experiment {other_exp}'] <= 0)
    
    exclusive_expression = experimental_data.loc[condition]
    exclusive_expression.to_csv(f"results/{exp}_exclusive_expression.csv", index=False)
    return exclusive_expression

# ===== MAIN FUNCTIONS =====

def process_with_controls(experimental_data_txt):
    """Process data for experiments with controls"""
    print("Normalize MaxQuant peptides.txt by Intensity and aggregate by Protein names\n")
    experiment_list = get_experiments_with_controls()
    experimental_data = convert_txt(experimental_data_txt)
    
    # Process data
    processed_data = remove_extra_columns_with_controls(experiment_list, experimental_data)
    
    # Run analyses
    print("\nRunning analyses...")
    get_differentially_expressed(experiment_list, processed_data)
    
    # Generate exclusivity reports for each experiment
    for i in range(len(experiment_list)):
        get_experiment_exclusive(experiment_list, processed_data, i)
    
    print("All analyses complete.")
    return processed_data

def process_without_controls(experimental_data_txt):
    """Process data for experiments without controls"""
    print("Filter MaxQuant data and aggregate rows by protein name")
    experiment_list = get_experiments_without_controls()
    experimental_data = convert_txt(experimental_data_txt)
    
    # Process data
    processed_data = remove_extra_columns_without_controls(experiment_list, experimental_data)
    
    # Run analyses
    print("\nRunning analyses...")
    get_differentially_expressed(experiment_list, processed_data)
    
    # Generate exclusivity reports for each experiment
    for i in range(len(experiment_list)):
        get_experiment_exclusive(experiment_list, processed_data, i)
    
    print("All analyses complete.")
    return processed_data

def main():
    """Main function to run the program"""
    parser = argparse.ArgumentParser(description='Process MaxQuant proteomics data.')
    parser.add_argument('input_file', help='Path to the peptides.txt file from MaxQuant')
    args = parser.parse_args()
    
    print("===== MaxQuant Proteomics Data Analysis =====")
    print("This program will process your proteomics data and generate analysis reports.")
    
    analysis_type = get_analysis_type()
    
    if analysis_type == 1:
        # Process with controls
        process_with_controls(args.input_file)
    else:
        # Process without controls
        process_without_controls(args.input_file)
    
    print("\nDone. Results saved in the 'results' directory.")

if __name__ == "__main__":
    main()
