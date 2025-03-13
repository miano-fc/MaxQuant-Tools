# Import Modules
import pandas as pd
import numpy as np

# Get a list of experimental conditions (2)
def get_experiments():
    print("\n")
    experiment_1 = input("Enter experiment 1: ")
    experiment_2 = input("Enter experiment 2: ")

    return [experiment_1, experiment_2]

# Get a list of proteins that are differentially expressed and calculate Log2FC by intensity
def select_log2fc():
    selection = input("\nHow should Log2 fold change be calculated for intensity?  (Experiment 1/2 or 2/1): ")
    return selection

def get_differentially_expressed(experiment_list, experimental_data):
    print("Generating a list of differentially-expressed proteins.\n") 
    
    differential_expression = experimental_data.loc[
        (experimental_data['Count ' + experiment_list[0]] != 0) &
        (experimental_data['Count ' + experiment_list[1]] != 0) &
        (experimental_data['Intensity Experiment ' + experiment_list[0]] > 0) &
        (experimental_data['Intensity Experiment ' + experiment_list[1]] > 0) 
        ]

    user_selection = select_log2fc()

    if user_selection == "1/2":
        differential_expression.loc[:, 'Log2FC 1/2'] = np.log2(differential_expression['Intensity Experiment ' + experiment_list[0]] / differential_expression['Intensity Experiment ' + experiment_list[1]])

    if user_selection == "2/1":
        differential_expression.loc[:, 'Log2FC 2/1'] = np.log2(differential_expression['Intensity Experiment ' + experiment_list[1]] / differential_expression['Intensity Experiment ' + experiment_list[0]])

    differential_expression.to_csv("results/differential_expression.csv", index = False)

    return differential_expression

# Get a lsit of proteins that are only expressed in experiment 1
def get_experiment_1(experiment_list, experimental_data):
    print("Generating a list of proteins expressed in experiment 1.")

    experiment_1_expression = experimental_data.loc[
        (experimental_data['Count ' + experiment_list[0]] > 0) & 
        (experimental_data['Count ' + experiment_list[1]] <= 0) &
        (experimental_data['Intensity Experiment ' + experiment_list[0]] > 0) &
        (experimental_data['Intensity Experiment ' + experiment_list[1]] <= 0)
        ]

    experiment_1_expression.to_csv("results/experiment_1_expression.csv", index = False)
    
    return experiment_1_expression

# Get a list of proteins that are only expressed in experiment 2
def get_experiment_2(experiment_list, experimental_data):
    print("Generating a list of proteins expressed in experiment 2.")

    experiment_2_expression = experimental_data.loc[
        (experimental_data['Count ' + experiment_list[1]] > 0) &
        (experimental_data['Count ' + experiment_list[0]] <= 0) &
        (experimental_data['Intensity Experiment ' + experiment_list[1]] > 0) &
        (experimental_data['Intensity Experiment ' + experiment_list[0]] <= 0)
        ]

    experiment_2_expression.to_csv("results/experiment_2_expression.csv", index = False)
    
    return experiment_2_expression

# Define a main function
def main(experimental_data):
    experiment_list = get_experiments()

    get_differentially_expressed(experiment_list, experimental_data)
    get_experiment_1(experiment_list, experimental_data)
    get_experiment_2(experiment_list, experimental_data)

    return get_differentially_expressed, get_experiment_1, get_experiment_2
