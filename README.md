# MaxQuant Tools
### What's ✨new✨
Modified the program to allow users to specify how many experiments and controls they have

### What I'm working on 🛠️
- A GUI via tkinter
- Streamlining the application/reducing overlapping code

---

## Purpose
To compare results of two [MaxQuant](https://www.maxquant.org) experiments and their respective controls from a run by normalizing data from the output peptides.txt file.

## Functionality
Users can selected data by PEP and MS/MS count thresholds, proteins marked by MaxQuant as potential contaminats are removed by default.

Data is normalized by count and intensity columns: (Experiment 1 - Control 1, Experiment 2 - Control 2).

## Output 
1. A .csv file containing aggregated, normalized data
2. A .csv file containing differentially expressed proteins + log2 fold change of intensity
3. .csv files containing proteins expressed in either Experiment 1 or Experiment 2

## How to use
In your Python environment, run the sort_data_master program:
<pre>
  python sort_data_master.py path/to/peptides.txt
</pre>

### Contributors
[Michael Miano](mailto:Michael.Miano@fccc.edu)
