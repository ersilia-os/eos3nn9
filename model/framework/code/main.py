# imports
import os
import csv
import sys
import pandas as pd
import subprocess
import base64
import pickle
import numpy as np
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt

np.bool = np.bool_
np.object = np.object_

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
root = os.path.dirname(os.path.abspath(__file__))

# my model
def my_model(smiles_list):
    return [MolWt(Chem.MolFromSmiles(smi)) for smi in smiles_list]


# Molecular descriptor calculator option
def desc_calc():
    # Performs the descriptor calculation
    bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/MACCSFingerprinter.xml -dir ./ -file descriptors_output.csv"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    # Read the original input file to get the order
    with open("molecule.smi", "r") as f:
        input_order = [line.split()[1] for line in f]

    # Read the output CSV
    descriptors_df = pd.read_csv("descriptors_output.csv")

    # Sort the output based on the input order
    sorted_descriptors_df = descriptors_df.set_index('Name').reindex(input_order).reset_index()
    
    # Save the sorted output to a new file
    sorted_descriptors_df.to_csv("sorted_descriptors_output.csv", index=False)

# File download option
def filedownload(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href

# Model building section
def build_model(input_data, load_data): # added load_data as param
    load_model = pickle.load(open('Mpro_model.pkl', 'rb'))
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    prediction_output = pd.Series(prediction, name='pIC50')
    molecule_name = pd.Series(load_data[1], name='molecule_name')
    df = pd.concat([molecule_name, prediction_output], axis=1)
    print(df) # display output
 


# read SMILES from .csv file, assuming one column with header
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]
#    ids = [f"smi_{x} for x in range(len(smiles_list)] # added
#    df = pd.DataFrame({"id":ids, "smiles":smiles}) # added

# write molecule.smi with tabs instead of spaces for Padel functionality
output_file = "molecule.smi"
with open(output_file, "w") as exportFile:
    for line in smiles_list:
        new_line = line.replace(' ', '\t')
        exportFile.write(new_line + '\n')

# define load_data for error 07.11.24
load_data = pd.read_table(output_file, sep='\t', header=None)

# run model
desc_calc()
desc = pd.read_csv('sorted_descriptors_output.csv') # read new sorted output instead of randomly sorted descriptors_output
Xlist = list(pd.read_csv('lists_of_descriptor.csv').columns)
desc_subset = desc[Xlist]
build_model(desc_subset, load_data)
outputs = my_model(smiles_list)

#check input and output have the same length
input_len = len(smiles_list)
output_len = len(outputs)
assert input_len == output_len

# write output in a .csv file
with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow(["value"])  # header
    for o in outputs:
        writer.writerow([o])
