#!/usr/bin/env python
# coding: utf-8

import os
import glob
import csv
import json
import math
import sys
import shutil
import numpy as np
import pandas as pd
import re
import periodictable
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import warnings

""" Turning of this warning: DeprecationWarning: DataFrameGroupBy.apply operated on the grouping columns. 
    This behavior is deprecated, and in a future version of pandas the grouping columns will be excluded from the operation.
    Either pass `include_groups=False` to exclude the groupings or explicitly select the grouping columns after groupby to silence this warning.
    final_df = modified_df.groupby('Dump File').apply(process_molecule_target_2).reset_index(drop=True) """
warnings.filterwarnings("ignore", category=DeprecationWarning)

# Planck's constant
h_eV = 4.1357e-15     # Planck's constant (eV·s)
hbar = 1.05457266e-34 # planck's constant/2pi (J.s)

# Vibrational frequency and bond length data from Tersoff potentials
vibrational_data = {
    ('Be', 'H'): {'Nu_e_it': 8.61964e13, 'r0_it': 1.3380, 'Nu_e_tt': 1.44275e13, 'r0_tt': 2.03500, 'Nu_e_ii': 4.14678e14, 'r0_ii': 0.74144}, # Carolina's potential
    ('Be', 'D'): {'Nu_e_it': 6.39396e13, 'r0_it': 1.3380, 'Nu_e_tt': 1.44275e13, 'r0_tt': 2.03500, 'Nu_e_ii': 2.93221e14, 'r0_ii': 0.74144}, # Carolina's potential
    ('Be', 'T'): {'Nu_e_it': 5.45722e13, 'r0_it': 1.3380, 'Nu_e_tt': 1.44275e13, 'r0_tt': 2.03500, 'Nu_e_ii': 2.39414e14, 'r0_ii': 0.74144}, # Carolina's potential
    ('B', 'W') : {'Nu_e_it': 3.20960e13, 'r0_it': 1.4707, 'Nu_e_tt': 3.00500e13, 'r0_tt': 1.56530, 'Nu_e_ii': 6.02760e12, 'r0_ii': 2.06120}, # Antoine's potential
    ('B', 'Ar'): {'Nu_e_it': 0.00000000, 'r0_it': 0.0000, 'Nu_e_tt': 3.00500e13, 'r0_tt': 1.56530, 'Nu_e_ii': 0.00000000, 'r0_ii': 0.00000}, # Antoine's potential
    ('W', 'H') : {'Nu_e_it': 5.58144e13, 'r0_it': 1.7270, 'Nu_e_tt': 7.43982e12, 'r0_tt': 2.34095, 'Nu_e_ii': 4.14678e14, 'r0_ii': 0.74144}, # Juslin's potential
    ('W', 'D') : {'Nu_e_it': 3.95741e13, 'r0_it': 1.7270, 'Nu_e_tt': 7.43982e12, 'r0_tt': 2.34095, 'Nu_e_ii': 2.93221e14, 'r0_ii': 0.74144}, # Juslin's potential
    ('W', 'T') : {'Nu_e_it': 3.24259e13, 'r0_it': 1.7270, 'Nu_e_tt': 7.43982e12, 'r0_tt': 2.34095, 'Nu_e_ii': 2.39414e14, 'r0_ii': 0.74144}, # Juslin's potential
}

# Take target and ion as input, capitalize for consistency
target = input("Enter the target element symbol (e.g., 'B'): ").strip().capitalize()
ion = input("Enter the ion element symbol (e.g., 'W'): ").strip().capitalize()

# Retrieve atomic masses
try:
    mass_target = getattr(periodictable, target).mass
    mass_ion = getattr(periodictable, ion).mass
    print(f"\n\033[32mMass of {target}: {mass_target} u\033[0m")
    print(f"\033[32mMass of {ion}: {mass_ion} u\n\033[0m")
except AttributeError:
    print("\033[31m\nOne or both element symbols are invalid. Please check your input. u\033[0m")
    sys.exit(1)

# Get vibrational data
key = (target, ion) if (target, ion) in vibrational_data else (ion, target)
data = vibrational_data.get(key, None)

if data:
    # Swap the values dynamically if the key order is reversed
    if key != (target, ion):
        Nu_e_tt, r0_tt = data['Nu_e_ii'], data['r0_ii']
        Nu_e_it, r0_it = data['Nu_e_it'], data['r0_it']
        Nu_e_ii, r0_ii = data['Nu_e_tt'], data['r0_tt']
    else:
        Nu_e_it, r0_it = data['Nu_e_it'], data['r0_it']
        Nu_e_tt, r0_tt = data['Nu_e_tt'], data['r0_tt']
        Nu_e_ii, r0_ii = data['Nu_e_ii'], data['r0_ii']

    print(f"\033[32mVibrational data for {ion}->{target}:\033[0m")
    print(f"\033[32m  Nu_e_it = {Nu_e_it:.5e} Hz, r0_it = {r0_it:.5e}\033[0m")
    print(f"\033[32m  Nu_e_tt = {Nu_e_tt:.5e} Hz, r0_tt = {r0_tt:.5e}\033[0m")
    print(f"\033[32m  Nu_e_ii = {Nu_e_ii:.5e} Hz, r0_ii = {r0_ii:.5e}\033[0m")
else:
    Nu_e_tt, r0_tt = 0, 0
    Nu_e_it, r0_it = 0, 0
    Nu_e_ii, r0_ii = 0, 0
    print(f"\033[32mNo vibrational data found for {ion}->{target}:\033[0m")
    print(f"\033[32m  Nu_e_it = {Nu_e_it:.1e} Hz, r0_it = {r0_it:.1e}\033[0m")
    print(f"\033[32m  Nu_e_tt = {Nu_e_tt:.1e} Hz, r0_tt = {r0_tt:.1e}\033[0m")
    print(f"\033[32m  Nu_e_ii = {Nu_e_ii:.1e} Hz, r0_ii = {r0_ii:.1e}\033[0m")


# Initialize an empty dictionary to store directories
dirs = {}

# Loop to allow user to input multiple directories
while True:
    # Take directory name input
    dir_name = input("\nEnter the directory name (e.g., '300K-0deg') or type 'done' to start: ").strip()
    
    # Break the loop if the user is done
    if dir_name.lower() == 'done':
        break
    
    # Take the corresponding path input
    dir_path = input(f"\nEnter the path for {dir_name}: ").strip()
    
    # Store the directory in the dictionary
    dirs[dir_name] = dir_path

# Print the directories to verify
print("\n\033[32mDirectories entered:\033[0m")
for dir_name, dir_path in dirs.items():
    print(f"\033[32m{dir_name}: {dir_path}\033[0m")

print("\033[32m\nRunning...\033[0m")

# Creating the final json file for sputtering yield
output_json = f'sputtering_yields_{ion}{target}.json'

""" Function to count the total and physically sputtered atoms from the "sputtered.data" file """
def count_sputtered_atoms(file_path):
    physical_count = total_count = 0
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip():  # Skip empty lines
                total_count += 1
                if line.split()[-1] == '0':  # Check for physical sputtering
                    physical_count += 1
    return physical_count, total_count

""" Function to get the number of impacts from the "event.csv" file """
def count_impacts(file_path):
    with open(file_path, 'r') as file:
        return sum(1 for line in file)

"""Compute the number of sputtered atoms in a simulation"""
def ingress_egress(fname):
    # Skip empty files
    if os.path.getsize(fname) == 0:
        print(f"\033[31m\nSkipping empty file: {fname}\033[0m")
        return None

    try:
        rid = os.path.splitext(os.path.basename(fname))[0]
        df = pd.read_csv(fname, header=0).set_index("time")
        # Strip whitespace from column names
        df.columns = df.columns.str.strip()
        i1 = i2 = e1 = e2 = sputtered = 0
        # "event" is boolean True if at any point count is 1
        if event := df["c1"].any() or df["c2"].any():
            # compute ingress/egress
            diff1 = df["c1"].diff()
            diff2 = df["c2"].diff()
            # count ingress/egress
            i1 = diff1.loc[diff1 > 0].sum().astype("int32")
            i2 = diff2.loc[diff2 > 0].sum().astype("int32")
            e1 = -diff1.loc[diff1 < 0].sum().astype("int32")
            e2 = -diff2.loc[diff2 < 0].sum().astype("int32")
            # no. of sputtered as minimum between ingress/egress counts
            sputtered = min(i1, i2, e1, e2)
        # Get the value of the "seed"
        seed = df.iloc[-1, -1]
        return rid, i1, i2, e1, e2, sputtered, event, seed
    except Exception as e:
        print(f"\033[31m\nSkipping file {fname}: {e}\033[0m")
        return None

""" Function to run ingress_egress for all simulations and generate the "event.csv" file """
def run_ingress_egress(root_dir):
    # Number of cores to use for parallel processing
    number_of_cores = 30

    # Create a ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=number_of_cores) as executor:
        results = []
        for i in range(1, len([item for item in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, item))])+1):
            folder_path = os.path.join(root_dir, str(i))
            # Find all .csv files in the folder
            csv_files = glob.glob(os.path.join(folder_path, "*.csv"))

            # Check if there are any .csv files and process them
            for file_path in csv_files:
                if os.path.exists(folder_path) and os.path.isfile(file_path):
                    results.append(executor.submit(ingress_egress, file_path))

    # Create a list of DataFrames for each folder's results
    dfs = []
    for future in results:
        result = future.result()
        if result:  # Skip None results
            columns = ["rid", "i1", "i2", "e1", "e2", "sputtered", "event", "seed"]
            df = pd.DataFrame([result], columns=columns)
            dfs.append(df)

    # Concatenate the DataFrames into a single DataFrame if there are results
    if dfs:
        final_df = pd.concat(dfs, ignore_index=True)

        # Filter and write only rows with unique "seed" values
        unique_seed_df = final_df.drop_duplicates(subset="seed")

        # Output results to a CSV file
        output = os.path.join(root_dir, "event.csv")
        unique_seed_df.to_csv(output, index=False)
    else:
        print("\033[31m\nNo valid data to write to event.csv file.\033[0m")
        sys.exit(1)

""" Function to store all the properties of the sputtered target atoms, in "sputtered.data" file """
def generate_sputtered_data(root_dir):
        
    input_file = os.path.join(root_dir, "event.csv")
    output = os.path.join(root_dir, "sputtered.data")

    if os.path.exists(output):
        os.remove(output)

    # Read the CSV file
    df = pd.read_csv(input_file)

    # Filter rows where "sputtered" is non-zero
    non_zero_sputtered_df = df[df["sputtered"] != 0]

    # Extract the "rid" values corresponding to non-zero "sputtered"
    rid_values = non_zero_sputtered_df["rid"].tolist()

    # List to store rid_values where potential energy of sputtered target atom is not zero
    global rids
    rids = []
    
    # Iterate through the RID values
    for rid in rid_values:
        # Extract the numeric portion from the "rid" value
        rid_numeric = int(rid.split(".")[1])
    
        # Generate the folder path for the corresponding RID
        folder_path = os.path.join(root_dir, str(rid_numeric))
    
        # Check if the folder exists
        if os.path.exists(folder_path):
            # Create the text file name based on the RID
            txt_file_name = f"control.{rid_numeric}.txt"
            txt_file_path = os.path.join(folder_path, txt_file_name)
        
            # Check if the text file exists
            if os.path.exists(txt_file_path):
                # Read the text file and process lines containing the target atom with unique first column
                unique_first_columns = set()  # To store unique first columns
                output_lines = []  # To store lines that meet the requirements
            
                with open(txt_file_path, "r") as txt_file:
                    for line in txt_file:
                        if f' {target} ' in line:
                            columns = line.strip().split()
                        
                            # Check if the first column is unique
                            first_column = columns[0]
                            if first_column not in unique_first_columns:
                                unique_first_columns.add(first_column)
                                output_lines.append(f'rid: {rid_numeric} - {line}')
                            
                                # Check if potential energy is not zero
                                potential_energy = float(columns[-1])
                                if potential_energy != 0:
                                    rids.append(rid)
            
                with open(output, "a") as output_file:
                     output_file.writelines(output_lines)

            else:
                txt_file_name = f"test.{rid_numeric}.txt"
                txt_file_path = os.path.join(folder_path, txt_file_name)
        
                # Check if the text file exists
                if os.path.exists(txt_file_path):
                    # Read the text file and process lines containing the target atom with unique first column
                    unique_first_columns = set()  # To store unique first columns
                    output_lines = []  # To store lines that meet the requirements
            
                    with open(txt_file_path, "r") as txt_file:
                        for line in txt_file:
                            if f' {target} ' in line:
                                columns = line.strip().split()
                        
                                # Check if the first column is unique
                                first_column = columns[0]
                                if first_column not in unique_first_columns:
                                    unique_first_columns.add(first_column)
                                    output_lines.append(f'rid: {rid_numeric} - {line}')
                            
                                    # Check if potential energy is not zero
                                    potential_energy = float(columns[-1])
                                    if potential_energy != 0:
                                        rids.append(rid)
            
                    with open(output, "a") as output_file:
                         output_file.writelines(output_lines)
               
    return rids

""" Function to calculate distance between two atoms """
def calculate_distance(atom1, atom2):
    x1, y1, z1 = float(atom1['x']), float(atom1['y']), float(atom1['z'])
    x2, y2, z2 = float(atom2['x']), float(atom2['y']), float(atom2['z'])
    return math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)

""" Function to calculate velocity """
def calculate_velocity(x1, y1, z1, t1, x2, y2, z2, t2):
    # Velocity formula: v = (s2 - s1) / (t2 - t1)
    delta_x = x2 - x1
    delta_y = y2 - y1
    delta_z = z2 - z1
    delta_t = t2 - t1

    velocity_x = delta_x / delta_t
    velocity_y = delta_y / delta_t
    velocity_z = delta_z / delta_t

    overall_velocity = np.sqrt(velocity_x**2 + velocity_y**2 + velocity_z**2)

    return velocity_x, velocity_y, velocity_z, overall_velocity

""" Function to calculate kinetic energy """
def calculate_KE(m,v):
    ke = 1/2*((m*1.6605e-27)*(v*100)**2) # in J
    ke = ke * 6.242e+18 # in eV
    return ke

""" Function to process each group of ID within a dump file to clean the polyatomic_target.csv file """
def process_molecule_target_1(group):
    # Reset the index for each ID group
    group = group.reset_index(drop=True)

    # Check the number of unique ID_ion and ID_target values
    unique_id_ion_count = group[f'ID_{ion}'].nunique()
    unique_id_target_count = group[f'ID_{target}'].nunique()
    
    # If there are two unique ID_ion values, keep only the rows where the time step is repeated
    if unique_id_target_count == 1 and unique_id_ion_count == 1:
        selected_rows = group.head(2)

    else:
        # Find the repeated time steps within this dump file
        repeated_time_steps = group[group.duplicated(subset='Time Step')]['Time Step'].unique()

        # Keep only the rows with repeated time steps
        selected_rows = group[group['Time Step'].isin(repeated_time_steps)]
        
    return selected_rows

""" Define a function to process each group of Dump File to clean the polyatomic_target.csv file """
def process_molecule_target_2(group):
    # Reset the index for each Dump File group
    group = group.reset_index(drop=True)
    
    # Check the number of unique ID_ion and ID_target values
    unique_id_ion_count = group[f'ID_{ion}'].nunique()
    unique_id_target_count = group[f'ID_{target}'].nunique()
    
    # If there are only 2 rows, keep them
    if len(group) == 2:
        selected_rows = group
        
    elif len(group) > 2 and unique_id_ion_count == 1:
        # Identify the maximum repeated time steps
        time_step_counts = group['Time Step'].value_counts()
        max_repeated_count = time_step_counts.max()
        # Number of rows to keep
        a = max_repeated_count * 2

        if a < 7:
            # Keep rows with the time steps that have the maximum count
            selected_rows = group[group['Time Step'].isin(time_step_counts[time_step_counts == max_repeated_count].index)].head(a)
        elif a > 7:
            max_repeated_count = 3
            a = 6
            selected_rows = group[group['Time Step'].isin(time_step_counts[time_step_counts == max_repeated_count].index)].head(a)
            
    elif len(group) > 2 and unique_id_ion_count > 1:
        # Identify the maximum repeated time steps
        time_step_counts = group['Time Step'].value_counts()
        max_repeated_count = time_step_counts.max()
        # Number of rows to keep
        a = max_repeated_count * 2

        if a < 7:
            # Keep rows with the time steps that have the maximum count
            selected_rows = group[group['Time Step'].isin(time_step_counts[time_step_counts == max_repeated_count].index)].head(a)
        elif a > 7:
            max_repeated_count = 3
            a = 6
            selected_rows = group[group['Time Step'].isin(time_step_counts[time_step_counts == max_repeated_count].index)].head(a)

    return selected_rows

""" Function to find the spettered molecules and single ions leaving the surface """
def generate_molecule_data(root_dir):
    # Output files for different cases
    output_molecule_target = os.path.join(root_dir, "Target_molecules.csv")
    diatomic_output =  os.path.join(root_dir, "diatomic_target.csv")
    polyatomic_output = os.path.join(root_dir, "polyatomic_target.csv")
    output_molecule_ion = os.path.join(root_dir, "ion_molecules.csv")
    output_single_ion = os.path.join(root_dir, "ion_single.csv")

    # Write headers for CSV files
    with open(output_molecule_target, 'w', newline='') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerow([
            "Dump File",
            "Time Step",
            "Time",
            f"ID_{target}",
            f"ID_{ion}",
            "Bond length",
            f"PE {target}",
            f"PE {ion}",
            f"KE {target}",
            f"KE {ion}",
            f"X_{target}",
            f"Y_{target}",
            f"Z_{target}",
            f"X_{ion}",
            f"Y_{ion}",
            f"Z_{ion}",
            f"Vx_{target}",
            f"Vy_{target}",
            f"Vz_{target}",
            f"Vx_{ion}",
            f"Vy_{ion}",
            f"Vz_{ion}"
        ])

    with open(output_molecule_ion, 'w', newline='') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerow([
            "Dump File",
            "Time",
            f"ID_{ion}1",
            f"ID_{ion}2",
            "Bond length",
            f"PE {ion}1",
            f"PE {ion}2",
            f"KE {ion}1",
            f"KE {ion}2",
            f"X_{ion}1",
            f"Y_{ion}1",
            f"Z_{ion}1",
            f"X_{ion}2",
            f"Y_{ion}2",
            f"Z_{ion}2",
            f"Vx_{ion}1",
            f"Vy_{ion}1",
            f"Vz_{ion}1",
            f"Vx_{ion}2",
            f"Vy_{ion}2",
            f"Vz_{ion}2"
        ])

    with open(output_single_ion, 'w', newline='') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerow([
            "Dump File",
            "Time",
            "ID",
            "PE",
            "KE",
            "Vx",
            "Vy",
            "Vz"
        ])

    # Number of cores to use for parallel processing
    number_of_cores = 30

    # Create a ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=number_of_cores) as executor:

        # Process each RID value
        for rid in range(1, len([item for item in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, item))])+1):
            # Extract the numeric portion from the "rid" value
            #rid_numeric = int(rid.split(".")[1])
            rid_numeric = rid
            ids_target = []
            ids_ion    = []

            # Generate the folder path for the corresponding RID
            folder_path = os.path.join(root_dir, str(rid_numeric))

            # Check if the folder exists
            if os.path.exists(folder_path):
                # Process each dump file in the folder
                for file_name in os.listdir(folder_path):
                    if file_name.endswith('.txt'):
                        file_path = os.path.join(folder_path, file_name)

                        # Read the dump file
                        with open(file_path, 'r') as file:
                            lines = file.readlines()

                        # Initialize variables
                        timestep_data = []
                        current_timestep = {'timestep': None, 'time': None, 'atoms': []}
                        headers = None  # Initialize headers outside the loop

                        # Parse the dump file
                        iterator = iter(lines)
                        for line in iterator:
                            if line.startswith("ITEM: TIMESTEP"):
                                if current_timestep is not None:
                                    timestep_data.append(current_timestep)
                                #current_timestep = {'timestep': None, 'time': None, 'atoms': []}
                                current_timestep['timestep'] = int(next(iterator, '').strip())  # Get the timestep from the next line
                                
                            elif line.startswith("ITEM: TIME"):
                                current_timestep = {'timestep': None, 'time': None, 'atoms': []}
                                try:
                                    time_line = next(iterator, '')
                                    if time_line is not None:
                                        current_timestep['time'] = float(time_line.strip())  # Get the time from the next line
                                    else:
                                        current_timestep['time'] = None
                                except ValueError:
                                    current_timestep['time'] = None

                            elif line.startswith("ITEM: ATOMS"):
                                headers = line.split()[2:]
                            elif headers is not None and not line.startswith("ITEM"):
                                atom_data = line.split()

                                # Ensure that the atom_data has the same length as headers
                                if len(atom_data) == len(headers):
                                    atom = {headers[i]: float(atom_data[i]) if i >= 3 and atom_data[i].replace('.', '').isdigit() else atom_data[i] for i in range(len(headers))}
                                    current_timestep['atoms'].append(atom)

                        # Append the last timestep
                        #if current_timestep is not None:
                        #    timestep_data.append(current_timestep)

                        # Process Target_molecules.csv for target atoms
                        target_timesteps = [timestep for timestep in timestep_data if any(atom['element'] == f'{target}' for atom in timestep['atoms'])]

                        with open(output_molecule_target, 'a', newline='') as csv_file:
                            csv_writer = csv.writer(csv_file)

                            for timestep in target_timesteps:
                                target_atoms = [atom for atom in timestep['atoms'] if atom['element'] == f'{target}']

                                for target_atom in target_atoms:
                                    ion_atoms = [a for a in timestep['atoms'] if a['element'] == f'{ion}']

                                    if ion_atoms:
                                        for ion_atom in ion_atoms:
                                            distance = calculate_distance(target_atom, ion_atom)

                                            if target_atom['c_CMP_PE'] != 0 and distance <= 2 and float(ion_atom['c_CMP_PE']) > -2:
                                                ids_target.append(target_atom['id'])
                                                ids_ion.append(ion_atom['id'])
                                                csv_writer.writerow([
                                                    file_name,
                                                    str(timestep['timestep']),
                                                    str(timestep['time']),
                                                    str(target_atom['id']),
                                                    str(ion_atom['id']),
                                                    str(distance),
                                                    str(target_atom['c_CMP_PE']),
                                                    str(ion_atom['c_CMP_PE']),
                                                    str(target_atom['c_CMP_KE']),
                                                    str(ion_atom['c_CMP_KE']),
                                                    str(target_atom['x']),
                                                    str(target_atom['y']),
                                                    str(target_atom['z']),
                                                    str(ion_atom['x']),
                                                    str(ion_atom['y']),
                                                    str(ion_atom['z']),
                                                    str(target_atom['vx']),
                                                    str(target_atom['vy']),
                                                    str(target_atom['vz']),
                                                    str(ion_atom['vx']),
                                                    str(ion_atom['vy']),
                                                    str(ion_atom['vz'])
                                                ])
                                                
                                            else:
                                                for target_atom2 in target_atoms:

                                                    if target_atom2 != target_atom and target_atom2['id'] not in ids_target:
                                                        distance = calculate_distance(target_atom, target_atom2)

                                                        if target_atom['c_CMP_PE'] != 0 and distance <= 2:
                                                            ids_target.append(target_atom['id'])
                                                            csv_writer.writerow([
                                                                file_name,
                                                                str(timestep['timestep']),
                                                                str(timestep['time']),
                                                                str(target_atom['id']),
                                                                str(target + ': ' + target_atom2['id']),
                                                                str(distance),
                                                                str(target_atom['c_CMP_PE']),
                                                                str(target_atom2['c_CMP_PE']),
                                                                str(target_atom['c_CMP_KE']),
                                                                str(target_atom2['c_CMP_KE']),
                                                                str(target_atom['x']),
                                                                str(target_atom['y']),
                                                                str(target_atom['z']),
                                                                str(target_atom2['x']),
                                                                str(target_atom2['y']),
                                                                str(target_atom2['z']),
                                                                str(target_atom['vx']),
                                                                str(target_atom['vy']),
                                                                str(target_atom['vz']),
                                                                str(target_atom2['vx']),
                                                                str(target_atom2['vy']),
                                                                str(target_atom2['vz'])
                                                            ])

                                    else:
                                        for target_atom2 in target_atoms:

                                            if target_atom2 != target_atom and target_atom2['id'] not in ids_target:
                                                distance = calculate_distance(target_atom, target_atom2)

                                                if target_atom['c_CMP_PE'] != 0 and distance <= 3:
                                                    ids_target.append(target_atom['id'])
                                                    csv_writer.writerow([
                                                        file_name,
                                                        str(timestep['timestep']),
                                                        str(timestep['time']),
                                                        str(target_atom['id']),
                                                        str(target + ': ' + target_atom2['id']),
                                                        str(distance),
                                                        str(target_atom['c_CMP_PE']),
                                                        str(target_atom2['c_CMP_PE']),
                                                        str(target_atom['c_CMP_KE']),
                                                        str(target_atom2['c_CMP_KE']),
                                                        str(target_atom['x']),
                                                        str(target_atom['y']),
                                                        str(target_atom['z']),
                                                        str(target_atom2['x']),
                                                        str(target_atom2['y']),
                                                        str(target_atom2['z']),
                                                        str(target_atom['vx']),
                                                        str(target_atom['vy']),
                                                        str(target_atom['vz']),
                                                        str(target_atom2['vx']),
                                                        str(target_atom2['vy']),
                                                        str(target_atom2['vz'])
                                                    ])

                        # Process ion_molecules.csv for D2, H2, etc.
                        ion_timesteps = [timestep for timestep in timestep_data if any(atom['element'] == f'{ion}' for atom in timestep['atoms'])]

                        with open(output_molecule_ion, 'a', newline='') as csv_file:
                            csv_writer = csv.writer(csv_file)

                            for timestep in ion_timesteps:
                                ion_atoms = [atom for atom in timestep['atoms'] if atom['element'] == f'{ion}']

                                for ion_atom1 in ion_atoms:

                                    for ion_atom2 in ion_atoms:

                                        if ion_atom1 != ion_atom2 and ion_atom1['id'] not in ids_ion:
                                            distance = calculate_distance(ion_atom1, ion_atom2)

                                            if ion_atom1['c_CMP_PE'] != 0 and ion_atom2['c_CMP_PE'] == ion_atom1['c_CMP_PE'] and distance <= 1:
                                                ids_ion.append(ion_atom2['id'])
                                                csv_writer.writerow([
                                                    file_name,
                                                    str(timestep['time']),
                                                    str(ion_atom1['id']),
                                                    str(ion_atom2['id']),
                                                    str(distance),
                                                    str(ion_atom1['c_CMP_PE']),
                                                    str(ion_atom2['c_CMP_PE']),
                                                    str(ion_atom1['c_CMP_KE']),
                                                    str(ion_atom2['c_CMP_KE']),
                                                    str(ion_atom1['x']),
                                                    str(ion_atom1['y']),
                                                    str(ion_atom1['z']),
                                                    str(ion_atom2['x']),
                                                    str(ion_atom2['y']),
                                                    str(ion_atom2['z']),
                                                    str(ion_atom1['vx']),
                                                    str(ion_atom1['vy']),
                                                    str(ion_atom1['vz']),
                                                    str(ion_atom2['vx']),
                                                    str(ion_atom2['vy']),
                                                    str(ion_atom2['vz'])
                                                ])


                        # Process ion_single.csv
                        with open(output_single_ion, 'a', newline='') as csv_file:
                            csv_writer = csv.writer(csv_file)

                            for timestep in ion_timesteps:

                                for atom in timestep['atoms']:

                                    if atom['element'] == f'{ion}' and float(atom['c_CMP_PE']) == 0 and atom['id'] not in ids_ion:
                                        ids_ion.append(atom['id'])
                                        csv_writer.writerow([
                                            file_name,
                                            str(timestep['time']),
                                            str(atom['id']),
                                            str(atom['c_CMP_PE']),
                                            str(atom['c_CMP_KE']),
                                            str(atom['vx']),
                                            str(atom['vy']),
                                            str(atom['vz'])
                                        ])
                            continue
    
    # Load the generated Molecule.csv file into a DataFrame
    molecule_df = pd.read_csv(output_molecule_target)

    # Create two lists to store diatomic and polyatomic molecules
    diatomic_target = []
    polyatomic_target = []

    # Group by 'Dump File' and 'ID_W' to process each simulation separately
    grouped = molecule_df.groupby(['Dump File', 'Time Step'])

    for (dump_file, molecule_id), group in grouped:
        # Check if 'ID_ion' is '-' for rows where we only have target2 molecules
        if group[f'ID_{ion}'].astype(str).str.startswith(target + ': ').iloc[0]:
            # Count unique 'ID_W' for W2 molecules
            unique_target_count = group[f'ID_{target}'].nunique()
            if unique_target_count == 1:
                diatomic_target.append(group)
            else:
                polyatomic_target.append(group)
        else:
            # For cases where ID_ion is present, count unique target and ion atoms
            unique_atom_count = group[f'ID_{target}'].nunique() + group[f'ID_{ion}'].nunique()
            if unique_atom_count == 2:
                diatomic_target.append(group)
            else:
                polyatomic_target.append(group)

    # Concatenate the results back into DataFrames
    if diatomic_target:
        diatomic_df = pd.concat(diatomic_target)
        diatomic_df.to_csv(diatomic_output, index=False)
 

    if polyatomic_target:
        polyatomic_df = pd.concat(polyatomic_target)
        polyatomic_df.to_csv(polyatomic_output, index=False)

    # Apply the processing function to each group of Dump File and save the final DataFrame to CSV
    modified_df_target = molecule_df.groupby('Dump File').apply(process_molecule_target_1).reset_index(drop=True)
    final_df_target = modified_df_target.groupby('Dump File').apply(process_molecule_target_2).reset_index(drop=True)
    final_df_target.to_csv(output_molecule_target, index=False)
    
    # Read the modified CSV file into a DataFrame
    df_target = pd.read_csv(output_molecule_target)

    # Filter rows where there are two or more rows for the same target atom ID in a dump file
    filtered_df_target = df_target.groupby(['Dump File', f'ID_{target}']).filter(lambda x: len(x) >= 2)

    # Calculate center of mass for each unique combination of "Dump File" and f"ID_{target}"
    com_data = []
    for (dump_file,), group in filtered_df_target.groupby(['Dump File']):
        dump_file = dump_file.strip()  # Remove leading and trailing whitespaces

        # Get unique time steps
        unique_time_steps = group['Time Step'].unique()

        for time_step in unique_time_steps:
            time_step_group = group[group['Time Step'] == time_step]

            # Calculate center of mass for each unique target atom and ion ID in that time step
            unique_target_ids = time_step_group[f'ID_{target}'].unique()
            unique_ion_ids = time_step_group[f'ID_{ion}'].unique()

            if all(not str(id_).startswith(target + ': ') for id_ in unique_ion_ids):
                total_mass = mass_target * len(unique_target_ids) + mass_ion * len(unique_ion_ids)
                com_x = ((time_step_group[time_step_group[f'ID_{target}'].isin(unique_target_ids)][f'X_{target}'] * mass_target).min() + (time_step_group[time_step_group[f'ID_{ion}'].isin(unique_ion_ids)][f'X_{ion}'] * mass_ion).sum()) / total_mass
                com_y = ((time_step_group[time_step_group[f'ID_{target}'].isin(unique_target_ids)][f'Y_{target}'] * mass_target).min() + (time_step_group[time_step_group[f'ID_{ion}'].isin(unique_ion_ids)][f'Y_{ion}'] * mass_ion).sum()) / total_mass
                com_z = ((time_step_group[time_step_group[f'ID_{target}'].isin(unique_target_ids)][f'Z_{target}'] * mass_target).min() + (time_step_group[time_step_group[f'ID_{ion}'].isin(unique_ion_ids)][f'Z_{ion}'] * mass_ion).sum()) / total_mass
                time = time_step_group[time_step_group[f'ID_{target}'].isin(unique_target_ids)]['Time'].min()
                molecule = f'{target}{len(unique_target_ids)}{ion}{len(unique_ion_ids)}'
                if len(unique_target_ids) == 1:
                    id_target = unique_target_ids[0]

                    # Append data to the result list
                    com_data.append([dump_file, time_step, time, id_target, molecule, total_mass, com_x, com_y, com_z])

                if len(unique_target_ids) == 2:
                    id_target = f'{unique_target_ids[0]} and {unique_target_ids[1]}'

                    # Append data to the result list
                    com_data.append([dump_file, time_step, time, id_target, molecule, total_mass, com_x, com_y, com_z])
            else:
                total_mass = mass_target * len(unique_target_ids)
                com_x = ((time_step_group[time_step_group[f'ID_{target}'].isin(unique_target_ids)][f'X_{target}'] * mass_target).sum()) / total_mass
                com_y = ((time_step_group[time_step_group[f'ID_{target}'].isin(unique_target_ids)][f'Y_{target}'] * mass_target).sum()) / total_mass
                com_z = ((time_step_group[time_step_group[f'ID_{target}'].isin(unique_target_ids)][f'Z_{target}'] * mass_target).sum()) / total_mass
                time = time_step_group[time_step_group[f'ID_{target}'].isin(unique_target_ids)]['Time'].min()
                num_t = f'{target}2'
                numeric_part = re.search(r'\d+', unique_ion_ids[0]).group()
                id_target = f'{unique_target_ids[0]} and {numeric_part}'
                # Append data to the result list
                com_data.append([dump_file, time_step, time, id_target, num_t, total_mass, com_x, com_y, com_z])

    # Create a DataFrame for center of mass data
    com_df_target = pd.DataFrame(com_data, columns=['Dump File', 'Time Step', 'Time', f'ID_{target}', 'Molecule', 'Mass', 'Center of Mass X', 'Center of Mass Y', 'Center of Mass Z'])

    # Write the result to the output CSV file
    com_df_target.to_csv(output_molecule_target, index=False)

    # Calculate Velocity and generate Velocity.csv

    # Read the center of mass CSV file
    com_df_target = pd.read_csv(output_molecule_target)

    # Calculate velocity for each pair of rows with the same 'Dump File'
    velocity_data = []
    for dump_file, group in com_df_target.groupby('Dump File'):
        # Sort by 'Time Step' to ensure consecutive rows
        group = group.sort_values(by='Time Step')

        # Iterate over consecutive pairs of rows
        for i in range(len(group) - 1):
            row1 = group.iloc[i]
            row2 = group.iloc[i + 1]
            if row1[f'ID_{target}'] == row2[f'ID_{target}']:
                # Calculate velocity
                velocity_x, velocity_y, velocity_z, overall_velocity = calculate_velocity(
                    row1['Center of Mass X'], row1['Center of Mass Y'], row1['Center of Mass Z'], row1['Time'],
                    row2['Center of Mass X'], row2['Center of Mass Y'], row2['Center of Mass Z'], row2['Time'])

                KE = calculate_KE(row1['Mass'], overall_velocity)

                # Append data to the result list
                velocity_data.append([dump_file, row1[f'ID_{target}'], row1['Molecule'], row1['Mass'], overall_velocity, KE])

    # Create a DataFrame for velocity data
    velocity_df_target = pd.DataFrame(velocity_data, columns=['Dump File', f'ID_{target}', 'Molecule', 'Mass', 'Velocity_com', 'KE_com'])

    # Write the result to the output CSV file
    velocity_df_target.to_csv(output_molecule_target, index=False)

    if os.path.exists(diatomic_output):
        # Read the modified CSV file into a DataFrame
        df_target = pd.read_csv(diatomic_output)
        target_data = []
        for index, row in df_target.iterrows():
            dump_file = row['Dump File'].strip()  # Remove leading and trailing whitespaces

            # Calculate center of mass for the unique ID in that line
            target_id = [row[f'ID_{target}']]
            ion_id = [row[f'ID_{ion}']]
            
            if str(ion_id[0]).startswith(target + ': '):
                # Compute the total mass of the system
                total_mass = mass_target * 2

                # Calculate the x, y, and z-coordinate of the center of mass (COM) using the weighted average formula
                COMx = ((row[f'X_{target}'] * mass_target) + (row[f'X_{ion}'] * mass_target)) / total_mass
                COMy = ((row[f'Y_{target}'] * mass_target) + (row[f'Y_{ion}'] * mass_target)) / total_mass
                COMz = ((row[f'Z_{target}'] * mass_target) + (row[f'Z_{ion}'] * mass_target)) / total_mass

                # Calculate the x, y, and z-component of the velocity of the center of mass
                Vcom_x = ((row[f'Vx_{target}'] * mass_target) + (row[f'Vx_{ion}'] * mass_target)) / total_mass
                Vcom_y = ((row[f'Vy_{target}'] * mass_target) + (row[f'Vy_{ion}'] * mass_target)) / total_mass
                Vcom_z = ((row[f'Vz_{target}'] * mass_target) + (row[f'Vz_{ion}'] * mass_target)) / total_mass

                # Compute the magnitude of the center-of-mass velocity vector
                Vcom = np.sqrt(Vcom_x**2 + Vcom_y**2 + Vcom_z**2)

                # Total kinetic energy (translation + rotational + vibrational) from the dataset
                KE_tot = row[f'KE {target}'] + row[f'KE {ion}']

                # Compute the translational kinetic energy (convert to SI units and then to eV)
                KE_com = (1/2*((total_mass*1.6605e-27)*(Vcom*100)**2)) * 6.242e+18
 
                # Compute the position vector of the target and ion atom relative to the COM
                r1 = [row[f'X_{target}'] - COMx, row[f'Y_{target}'] - COMy, row[f'Z_{target}'] - COMz]
                r2 = [row[f'X_{ion}'] - COMx, row[f'Y_{ion}'] - COMy, row[f'Z_{ion}'] - COMz]

                # Compute the velocity of the target and ion relative to the center of mass
                v_rel1 = [row[f'Vx_{target}'] - Vcom_x, row[f'Vy_{target}'] - Vcom_y, row[f'Vz_{target}'] - Vcom_z]
                v_rel2 = [row[f'Vx_{ion}'] - Vcom_x, row[f'Vy_{ion}'] - Vcom_y, row[f'Vz_{ion}'] - Vcom_z]
            
                # Compute the moment of inertia tensor for the target and ion
                I1 = mass_target * (np.dot(r1, r1) * np.identity(3) - np.outer(r1, r1))
                I2 = mass_target * (np.dot(r2, r2) * np.identity(3) - np.outer(r2, r2))

                # Total moment of inertia tensor for the diatomic molecule
                I = I1 + I2

                # Compute the angular momentum vector for the target and ion
                L1 = mass_target * np.cross(r1, v_rel1)
                L2 = mass_target * np.cross(r2, v_rel2)

                # Total angular momentum of the system
                L = L1 + L2
            
                # Solve for angular velocity vector using the pseudo-inverse of the inertia tensor
                omega = np.dot(np.linalg.pinv(I), L)
            
                # Compute rotational kinetic energy in eV
                KE_rot = 0.5 * np.dot(omega, np.dot(I, omega)) * 1.036427e-4
                
                # Compute vibrational kinetic energy
                KE_vib = KE_tot - KE_rot - KE_com
                
                # Reduced mass of the diatomic molecule in kg
                mu = ((mass_target * mass_target) / (mass_target + mass_target)) * 1.6605e-27

                if r0_tt !=0:
                    # Compute the rotational quantum number J based on rotational kinetic energy
                    J_term = (2 * mu * (r0_tt * 1e-10)**2 * KE_rot) / (hbar**2 * 6.242e+18)
                    J_temp = (-1 + np.sqrt(1 + 4 * J_term)) / 2  # Solve for J(J+1)

                    # Round to nearest integer for J (ensuring non-negative values)
                    J = np.round(J_temp).astype(int)
                    if J >= 0:
                        J = J
                    else:
                        J = 0
                else:
                    J = '-'

                if Nu_e_tt != 0:
                    # Compute the vibrational quantum number n based on vibrational energy
                    n_temp = (2 * KE_vib) / (h_eV * Nu_e_tt) - 0.5

                    # Round to nearest integer for n (ensuring non-negative values)
                    n = np.round(n_temp).astype(int)
                    if n >= 0:
                        n = n
                    else:
                        n = 0
                else:
                    n = '-'

                molecule = f'{target}{len(target_id) + len(ion_id)}'
                id_molecule = f'{target_id} and {ion_id}'

                # Append data to the result list
                target_data.append([dump_file,
                                row['Time'],
                                id_molecule,
                                molecule,
                                KE_tot,
                                KE_com,
                                KE_rot,
                                KE_vib,
                                J,
                                n,
                                Vcom,
                                Vcom_x,
                                Vcom_y,
                                Vcom_z,
                                COMx,
                                COMy,
                                COMz
                                ])            

            else:
                total_mass = mass_target + mass_ion
                COMx = ((row[f'X_{target}'] * mass_target) + (row[f'X_{ion}'] * mass_ion)) / total_mass
                COMy = ((row[f'Y_{target}'] * mass_target) + (row[f'Y_{ion}'] * mass_ion)) / total_mass
                COMz = ((row[f'Z_{target}'] * mass_target) + (row[f'Z_{ion}'] * mass_ion)) / total_mass

                Vcom_x = ((row[f'Vx_{target}'] * mass_target) + (row[f'Vx_{ion}'] * mass_ion)) / total_mass
                Vcom_y = ((row[f'Vy_{target}'] * mass_target) + (row[f'Vy_{ion}'] * mass_ion)) / total_mass
                Vcom_z = ((row[f'Vz_{target}'] * mass_target) + (row[f'Vz_{ion}'] * mass_ion)) / total_mass

                Vcom = np.sqrt(Vcom_x**2 + Vcom_y**2 + Vcom_z**2)

                KE_tot = row[f'KE {target}'] + row[f'KE {ion}']

                KE_com = (1/2*((total_mass*1.6605e-27)*(Vcom*100)**2)) * 6.242e+18

                r1 = [row[f'X_{target}'] - COMx, row[f'Y_{target}'] - COMy, row[f'Z_{target}'] - COMz]
                r2 = [row[f'X_{ion}'] - COMx, row[f'Y_{ion}'] - COMy, row[f'Z_{ion}'] - COMz]

                v_rel1 = [row[f'Vx_{target}'] - Vcom_x, row[f'Vy_{target}'] - Vcom_y, row[f'Vz_{target}'] - Vcom_z]
                v_rel2 = [row[f'Vx_{ion}'] - Vcom_x, row[f'Vy_{ion}'] - Vcom_y, row[f'Vz_{ion}'] - Vcom_z]
            
                # Moment of Inertia Tensor (I)       
                I1 = mass_target * (np.dot(r1, r1) * np.identity(3) - np.outer(r1, r1))
                I2 = mass_ion * (np.dot(r2, r2) * np.identity(3) - np.outer(r2, r2))
                I = I1 + I2
                
                # Angular Momentum (L): L = m*r*v
                L1 = mass_target * np.cross(r1, v_rel1)
                L2 = mass_ion * np.cross(r2, v_rel2)
                L = L1 + L2
            
                # omega = I^-1 . L
                omega = np.dot(np.linalg.pinv(I), L)
            
                KE_rot = 0.5 * np.dot(omega, np.dot(I, omega)) * 1.036427e-4
                
                KE_vib = KE_tot - KE_rot - KE_com
                
                mu = ((mass_target * mass_ion) / (mass_target + mass_ion)) * 1.6605e-27

                if r0_it != 0:
                    J_term = (2 * mu * (r0_it * 1e-10)**2 * KE_rot) / (hbar**2 * 6.242e+18)
                    J_temp = (-1 + np.sqrt(1 + 4 * J_term)) / 2  # Solve for J(J+1)
                    J = np.round(J_temp).astype(int)
                    if J >= 0:
                        J = J
                    else:
                        J = 0
                else:
                    J = '-'

                if Nu_e_it != 0:
                    n_temp = (2 * KE_vib) / (h_eV * Nu_e_it) - 0.5
                    n = np.round(n_temp).astype(int)
                    if n >= 0:
                        n = n
                    else:
                        n = 0
                else:
                    n = '-'

                molecule = f'{target}{len(target_id)}{ion}{len(ion_id)}'
                id_molecule = f'{target_id} and {ion_id}'

                # Append data to the result list
                target_data.append([dump_file,
                                row['Time'],
                                id_molecule,
                                molecule,
                                KE_tot,
                                KE_com,
                                KE_rot,
                                KE_vib,
                                J,
                                n,
                                Vcom,
                                Vcom_x,
                                Vcom_y,
                                Vcom_z,
                                COMx,
                                COMy,
                                COMz
                                ])

        # Create a DataFrame for center of mass data
        com_df_target = pd.DataFrame(target_data, columns=['Dump File',
                                                    'Time',
                                                    'ID',
                                                    'Molecule',
                                                    'KE_tot',
                                                    'KE_com',
                                                    'KE_rot',
                                                    'KE_vib',
                                                    'Rot quantum #',
                                                    'Vib quantum #',
                                                    'Vcom', 
                                                    'Vcom_x',
                                                    'Vcom_y',
                                                    'Vcom_z',
                                                    'COMx',
                                                    'COMy',
                                                    'COMz'
                                                    ])

        # Write the result to the output CSV file
        com_df_target.to_csv(diatomic_output, index=False)

        if r0_it != 0:

            final_rovib_data =  os.path.join(root_dir, "final_rovib_data.csv")  # Replace with your desired output file path
            rovib_data = pd.read_csv(diatomic_output)

            # Columns to average
            columns_to_average = ["Rot quantum #", "Vib quantum #"]

            # Group by 'Dump File' and 'ID', and compute the mean for numeric columns
            grouped_rovib = rovib_data.groupby(['Dump File', 'ID'])

            # Compute the mean for the rotational and vibrational quantum numbers
            averaged = grouped_rovib[columns_to_average].mean()
            averaged = averaged.round(0)

            # Extract the last row of each group for all other columns
            last_rows = grouped_rovib.last()

            # Combine the averaged and last row data
            rovib = last_rows.copy()
            rovib[columns_to_average] = averaged

            # Save the result to a new CSV file
            rovib.reset_index().to_csv(final_rovib_data, index=False)

    # Calculate center of mass for each unique combination of "Dump File" and f"ID_{ion}1"
    df_ion    = pd.read_csv(output_molecule_ion)
    ion_data = []
    for index, row in df_ion.iterrows():
        dump_file = row['Dump File'].strip()  # Remove leading and trailing whitespaces

        # Calculate center of mass for the unique ID in that line
        ion_id1 = [row[f'ID_{ion}1']]
        ion_id2 = [row[f'ID_{ion}2']]

        total_mass = mass_ion * (len(ion_id1) + len(ion_id2))
        COMx = ((row[f'X_{ion}1'] * mass_ion) + (row[f'X_{ion}2'] * mass_ion)) / total_mass
        COMy = ((row[f'Y_{ion}1'] * mass_ion) + (row[f'Y_{ion}2'] * mass_ion)) / total_mass
        COMz = ((row[f'Z_{ion}1'] * mass_ion) + (row[f'Z_{ion}2'] * mass_ion)) / total_mass

        Vcom_x = ((row[f'Vx_{ion}1'] * mass_ion) + (row[f'Vx_{ion}2'] * mass_ion)) / total_mass
        Vcom_y = ((row[f'Vy_{ion}1'] * mass_ion) + (row[f'Vy_{ion}2'] * mass_ion)) / total_mass
        Vcom_z = ((row[f'Vz_{ion}1'] * mass_ion) + (row[f'Vz_{ion}2'] * mass_ion)) / total_mass

        Vcom = np.sqrt(Vcom_x**2 + Vcom_y**2 + Vcom_z**2)

        KE_tot = row[f'KE {ion}1'] + row[f'KE {ion}2']

        KE_com = (1/2*((total_mass*1.6605e-27)*(Vcom*100)**2)) * 6.242e+18

        r1 = [row[f'X_{ion}1'] - COMx, row[f'Y_{ion}1'] - COMy, row[f'Z_{ion}1'] - COMz]
        r2 = [row[f'X_{ion}2'] - COMx, row[f'Y_{ion}2'] - COMy, row[f'Z_{ion}2'] - COMz]

        v_rel1 = [row[f'Vx_{ion}1'] - Vcom_x, row[f'Vy_{ion}1'] - Vcom_y, row[f'Vz_{ion}1'] - Vcom_z]
        v_rel2 = [row[f'Vx_{ion}2'] - Vcom_x, row[f'Vy_{ion}2'] - Vcom_y, row[f'Vz_{ion}2'] - Vcom_z]
    
        I1 = mass_ion * (np.dot(r1, r1) * np.identity(3) - np.outer(r1, r1))
        I2 = mass_ion * (np.dot(r2, r2) * np.identity(3) - np.outer(r2, r2))
        I = I1 + I2

        L1 = mass_ion * np.cross(r1, v_rel1)
        L2 = mass_ion * np.cross(r2, v_rel2)
        L = L1 + L2
    
        omega = np.dot(np.linalg.pinv(I), L)
    
        KE_rot = 0.5 * np.dot(omega, np.dot(I, omega)) * 1.036427e-4
        
        KE_vib = KE_tot - KE_rot - KE_com
        
        mu = ((mass_ion * mass_ion) / (mass_ion + mass_ion)) * 1.6605e-27

        J_term = (2 * mu * (float(row['Bond length']) * 1e-10)**2 * KE_rot) / (hbar**2 * 6.242e+18)
        l = (-1 + np.sqrt(1 + 4 * J_term)) / 2  # Solve for J(J+1)
        J = np.round(l).astype(int)
        if J >= 0:
            J = J
        else:
            J = 0

        k = (2 * KE_vib) / (h_eV * Nu_e_ii) - 0.5
        n = np.round(k).astype(int)
        if n >= 0:
            n = n
        else:
            n = 0

        molecule = f'{ion}{len(ion_id1) + len(ion_id2)}'
        id_ion = f'{ion_id1} and {ion_id2}'

        # Append data to the result list
        ion_data.append([dump_file,
                         row['Time'],
                         id_ion,
                         molecule,
                         KE_tot,
                         KE_com,
                         KE_rot,
                         KE_vib,
                         J,
                         n,
                         Vcom,
                         Vcom_x,
                         Vcom_y,
                         Vcom_z,
                         COMx,
                         COMy,
                         COMz
                        ])

    # Create a DataFrame for center of mass data
    com_df_ion = pd.DataFrame(ion_data, columns=['Dump File',
                                                 'Time',
                                                 f'ID_{ion}',
                                                 'Molecule',
                                                 'KE_tot',
                                                 'KE_com',
                                                 'KE_rot',
                                                 'KE_vib',
                                                 'Rot quantum #',
                                                 'Vib quantum #',
                                                 'Vcom', 
                                                 'Vcom_x',
                                                 'Vcom_y',
                                                 'Vcom_z',
                                                 'COMx',
                                                 'COMy',
                                                 'COMz'
                                                ])

    # Write the result to the output CSV file
    com_df_ion.to_csv(output_molecule_ion, index=False)

""" Function to calculate the type and number of sputtered species"""
def sputtered_species(name, root_dir):
    
    values_to_count = {
        f'{target}1{ion}1': f'{target}{ion}', 
        f'{target}1{ion}2': f'{target}{ion}2', 
        f'{target}1{ion}3': f'{target}{ion}3',
        f'{target}1{ion}4': f'{target}{ion}4',  
        f'{target}2'      : f'{target}2', 
        f'{target}2{ion}1': f'{target}2{ion}',
        f'{target}2{ion}2': f'{target}2{ion}2',
        f'{target}2{ion}3': f'{target}2{ion}3',
        f'{target}2{ion}4': f'{target}2{ion}4',
        f'{target}2{ion}5': f'{target}2{ion}5',
    }
    
    # Initialize results list
    results1 = []
    results2 = []
    energiez = []
    
    if os.path.isdir(root_dir):
        for E_dir in os.listdir(root_dir):
            E_path = os.path.join(root_dir, E_dir)
            if os.path.isdir(E_path) and E_dir.startswith(f'{ion}'):
                energy = int(E_dir[len(ion):])  # Extract energy value as integer
                energiez.append(energy)
                energiez.sort()
                folder_path = os.path.join(E_path, 'Target_molecules.csv')
                sputtered_data_path = os.path.join(E_path, 'sputtered.data')
                if os.path.exists(folder_path):
                    df = pd.read_csv(folder_path)
                    
                    for value, label in values_to_count.items():
                        count = df[df['Molecule'] == value].shape[0]
                        results1.append({'Energy [eV]' : energy, 'Value': label, 'Count': count})
                else:
                    #print(f"\n{name}: 'Velocity.csv' not found.")
                    for value, label in values_to_count.items():
                        count = 0
                        results1.append({'Energy [eV]' : energy, 'Value': label, 'Count': count})
                
                if os.path.exists(sputtered_data_path):
                    with open(sputtered_data_path, 'r') as file:
                        lines = file.readlines()
                        total_lines = len(lines)
                        results2.append({'Energy [eV]' : energy, 'Total Sputtered Species': total_lines})
                else:
                    #print(f"\n{name}: 'sputtered.data' not found.")
                    total_lines = 0
                    results2.append({'Energy [eV]' : energy, 'Total Sputtered Species': total_lines})
    
    # Create a DataFrame from the results list
    results1_df = pd.DataFrame(results1)
    results2_df = pd.DataFrame(results2)
    
    # Calculate sums for each value across all folders
    sums1_df = results1_df.groupby(['Energy [eV]', 'Value'])['Count'].sum().reset_index()
    sums2_df = results2_df.groupby(['Energy [eV]'])['Total Sputtered Species'].sum().reset_index()

    # Pivot the DataFrame to have molecule types as columns
    pivot1_df = sums1_df.pivot(index='Energy [eV]', columns='Value', values='Count').reset_index()
    pivot2_df = pd.DataFrame(sums2_df['Total Sputtered Species'].values, index=sums2_df['Energy [eV]'], columns=['Total Sputtered Species']).reset_index()
    
    # Calculate the 'target' column dynamically
    pivot1_df[f'{target}'] = pivot2_df['Total Sputtered Species'] - pivot1_df[[f'{target}{ion}', f'{target}{ion}2', f'{target}{ion}3', f'{target}2', f'{target}2', f'{target}2{ion}', f'{target}2{ion}']].sum(axis=1)
    
    # Sort the DataFrame by 'Folder'
    pivot1_df['Energy [eV]'] = pd.Categorical(pivot1_df['Energy [eV]'], categories=energiez, ordered=True)
    pivot1_df = pivot1_df.sort_values(by='Energy [eV]')

    # Save the results to a CSV file
    pivot1_df.to_csv(f'{ion}{target}_{name}_sputtered_species.csv', index=False)

################################################################################################################

# Function to process each directory concurrently
def process_directory(name, dire, ion):
    energies = set()
    results = {}

    # Your existing logic for directory analysis
    if os.path.isdir(dire):
        for energy_dir in os.listdir(dire):
            energy_path = os.path.join(dire, energy_dir)
            if os.path.isdir(energy_path) and energy_dir.startswith(f'{ion}'):
                energy = int(energy_dir[len(ion):])  # Extract energy value as integer
                energies.add(energy)
                sputtered_file = os.path.join(energy_path, 'sputtered.data')
                event_file = os.path.join(energy_path, 'event.csv')

                # If necessary files don't exist, generate them
                if not os.path.exists(event_file):
                    run_ingress_egress(energy_path)
                    generate_sputtered_data(energy_path)
                    generate_molecule_data(energy_path)

                if os.path.exists(event_file):
                    if os.path.exists(sputtered_file):
                        physical_count, total_count = count_sputtered_atoms(sputtered_file)
                    else:
                        physical_count = 0
                        total_count = 0

                    total_events = count_impacts(event_file) - 1

                    if total_events > 0:
                        physical_yield = physical_count / total_events
                        physical_error = math.sqrt(abs((physical_yield * (1 - physical_yield)) / total_events))
                        total_yield = total_count / total_events
                        total_error = math.sqrt(abs(total_yield * (1 - total_yield)) / total_events)
                        results[energy] = {
                            'physical yield': physical_yield,
                            'physical error': physical_error,
                            'total yield': total_yield,
                            'total error': total_error
                        }
                    else:
                        results[energy] = {
                            'physical yield': 'N/A (event.csv is empty)',
                            'physical error': 'N/A (event.csv is empty)',
                            'total yield': 'N/A (event.csv is empty)',
                            'total error': 'N/A (event.csv is empty)'
                        }
                else:
                    results[energy] = {
                        'physical yield': 'N/A (event.csv is missing)',
                        'physical error': 'N/A (event.csv is missing)',
                        'total yield': 'N/A (event.csv is missing)',
                        'total error': 'N/A (event.csv is missing)'
                    }
    
    sputtered_species(name, dire)

    return name, results, energies

# Main function to execute analysis in parallel
def main():
    # Prepare data structures for storing results
    energies = set()
    yields = {temp: {} for temp in dirs.keys()}

    # Use ThreadPoolExecutor to parallelize the directory processing
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        future_to_dir = {executor.submit(process_directory, name, dire, ion): name for name, dire in dirs.items()}

        for future in as_completed(future_to_dir):
            name, results, dir_energies = future.result()
            energies.update(dir_energies)
            yields[name] = results

    # Sort energies
    sorted_energies = sorted(energies)

    # Convert the yields dictionary to a format suitable for JSON
    data_for_json = {'Yield': {}}

    for dir_name in dirs.keys():
        data_for_json['Yield'][dir_name] = []
        for energy in sorted_energies:
            if energy in yields[dir_name]:
                physical_yield_value = yields[dir_name][energy]['physical yield']
                physical_error_value = yields[dir_name][energy]['physical error']
                total_yield_value = yields[dir_name][energy]['total yield']
                total_error_value = yields[dir_name][energy]['total error']

                if isinstance(physical_yield_value, float) and isinstance(physical_error_value, float):
                    data_for_json['Yield'][dir_name].append({
                        'Energy': energy,
                        'Physical Yield': f'{physical_yield_value:.6f}',
                        'Physical Error': f'{physical_error_value:.6f}',
                        'Total Yield': f'{total_yield_value:.6f}',
                        'Total Error': f'{total_error_value:.6f}'
                    })
                else:
                    data_for_json['Yield'][dir_name].append({
                        'Energy': energy,
                        'Physical Yield': physical_yield_value,
                        'Total Yield': total_yield_value
                    })

    # Write the results to the output JSON file
    with open(output_json, 'w') as jsonfile:
        json.dump(data_for_json, jsonfile, indent=4)

    print(f"\033[32m\nResults saved to {output_json}\033[0m")

if __name__ == "__main__":
    main()