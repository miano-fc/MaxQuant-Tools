import pandas as pd
import numpy as np

def get_experiments():
    experiments = []
    print("\nEnter experimental conditions (type 'done' when finished):")
    while True:
        experiment = input(f"Enter experiment {len(experiments) + 1}: ")
        if experiment.lower() == 'done':
            break
        experiments.append(experiment)
    return experiments

def select_log2fc(experiment_list):
    print("\nSelect how Log2 fold change should be calculated:")
    print("Options:")
    for i in range(len(experiment_list)):
        for j in range(i + 1, len(experiment_list)):
            print(f"  {i + 1}/{j + 1} ( {experiment_list[i]} / {experiment_list[j]} )")
            print(f"  {j + 1}/{i + 1} ( {experiment_list[j]} / {experiment_list[i]} )")
    
    selection = input("Enter your choice (e.g., 1/2 or 2/1): ")
    return selection

def get_differentially_expressed(experiment_list, experimental_data):
    print("Generating a list of differentially expressed proteins.\n")
    
    condition = (experimental_data[[f'Count {exp}' for exp in experiment_list]] != 0).all(axis=1) & \
                (experimental_data[[f'Intensity Experiment {exp}' for exp in experiment_list]] > 0).all(axis=1)
    differential_expression = experimental_data.loc[condition]
    
    user_selection = select_log2fc(experiment_list)
    exp1_idx, exp2_idx = map(int, user_selection.split('/'))
    exp1, exp2 = experiment_list[exp1_idx - 1], experiment_list[exp2_idx - 1]
    
    differential_expression[f'Log2FC {exp1}/{exp2}'] = np.log2(
        differential_expression[f'Intensity Experiment {exp1}'] /
        differential_expression[f'Intensity Experiment {exp2}']
    )
    
    differential_expression.to_csv("results/differential_expression.csv", index=False)
    return differential_expression

def get_experiment_exclusive(experiment_list, experimental_data, exp_index):
    exp = experiment_list[exp_index]
    print(f"Generating a list of proteins expressed exclusively in {exp}.")
    
    condition = (experimental_data[f'Count {exp}'] > 0) & \
                (experimental_data[f'Intensity Experiment {exp}'] > 0)
    
    for i, other_exp in enumerate(experiment_list):
        if i != exp_index:
            condition &= (experimental_data[f'Count {other_exp}'] <= 0) & \
                         (experimental_data[f'Intensity Experiment {other_exp}'] <= 0)
    
    exclusive_expression = experimental_data.loc[condition]
    exclusive_expression.to_csv(f"results/{exp}_exclusive_expression.csv", index=False)
    return exclusive_expression

def main(experimental_data):
    experiment_list = get_experiments()
    
    get_differentially_expressed(experiment_list, experimental_data)
    for i in range(len(experiment_list)):
        get_experiment_exclusive(experiment_list, experimental_data, i)

    print("Done.")
    
    return get_differentially_expressed, get_experiment_exclusive
