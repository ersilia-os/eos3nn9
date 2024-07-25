# imports
import os
import csv
import sys
import pandas as pd
import subprocess
import base64
import pickle
import numpy as np
import tempfile
import shutil

root = os.path.dirname(os.path.abspath(__file__))

padel_folder = os.path.abspath(os.path.join(root, "..", "PaDEL-Descriptor"))
checkpoints_folder = os.path.abspath(os.path.join(root, "..", "..", "checkpoints"))
tmp_folder = tempfile.mkdtemp(prefix="ersilia-")

np.bool = np.bool_
np.object = np.object_

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]


def desc_calc():
    # Performs the descriptor calculation with PADEL
    bashCommand = f"java -Xms2G -Xmx2G -Djava.awt.headless=true -jar {padel_folder}/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes {padel_folder}/MACCSFingerprinter.xml -dir {tmp_folder} -file {tmp_folder}/descriptors_output.csv"
    print(bashCommand)
    subprocess.Popen(bashCommand, shell=True).wait()
    #process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    #output, error = process.communicate()

    # Read the output CSV
    descriptors_df = pd.read_csv( f"{tmp_folder}/descriptors_output.csv" )
    # Sort the output based on the input order
    sorted_descriptors_df = descriptors_df.sort_values( by=descriptors_df.columns[0] )
    sorted_descriptors_df.to_csv( f"{tmp_folder}/sorted_descriptors_output.csv", index=False )

# Model building section
def build_model( input_data ):
    load_model = pickle.load(open(f'{checkpoints_folder}/Mpro_model.pkl', 'rb'))
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    df = pd.Series(prediction, name='pIC50')
    return( df ) 

# read SMILES from .csv file, assuming one column with header
ids = [] # declare ids list
with open( input_file, "r" ) as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]

# write molecule.smi for PADEL functionality
padel_file = f"{tmp_folder}/molecule.smi"
with open( padel_file, "w" ) as exportFile:
   for line in smiles_list:
        new_line = line.replace(' ', '\t') # not necessary for single column file
        exportFile.write(new_line + '\n')

# run model
desc_calc()
desc = pd.read_csv(f'{tmp_folder}/sorted_descriptors_output.csv')
Xlist = list(pd.read_csv(f'{checkpoints_folder}/lists_of_descriptor.csv').columns)
desc_subset = desc[Xlist]
outputs = build_model( desc_subset )
#check input and output have the same length
input_len = len(smiles_list)
output_len = len(outputs)
assert input_len == output_len

# write output in a .csv file
with open( output_file.strip(), "w") as f:
    writer = csv.writer(f)
    writer.writerow(["value"])  # header
    for o in outputs:
        writer.writerow([o])

shutil.rmtree(tmp_folder)