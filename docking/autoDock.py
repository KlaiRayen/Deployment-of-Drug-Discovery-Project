import subprocess
import re

def find_ligand_coordinates(chimerax_path, pdb_file, ligand_identifier, output_txt_file):
    chimerax_command = (
        f'open {pdb_file}; '
        f'measure center /a:{ligand_identifier}; '
        f'close all; exit'
    )
    command = [chimerax_path, '--nogui', '--cmd', chimerax_command]
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print(result.stdout)
        # Extract coordinates from the output
        match = re.search(r"Center of mass of \d+ atoms = \(([-\d.]+), ([-\d.]+), ([-\d.]+)\)", result.stdout)
        if match:
            center_x, center_y, center_z = match.groups()
            with open(output_txt_file, 'w') as file:
                file.write(f"center_x = {center_x}\n")
                file.write(f"center_y = {center_y}\n")
                file.write(f"center_z = {center_z}\n")
        else:
            print("Could not find the center of mass in the output.")
    except subprocess.CalledProcessError as e:
        print(f"Command failed with error: {e.stderr}")

def run_vina(config_file, output_file):
    command = [
        'C:/Users/klair/Desktop/DockingJob/Vina/vina.exe',  # Full path to vina executable
        '--config', config_file,
        '--out', output_file
    ]
    subprocess.run(command, check=True)

if __name__ == "__main__":
    # Define paths and parameters
    pdb_file = 'C://Users//klair//Desktop//AutomatedDocking//Automated-Docking//docking//1stp.pdb'
    chimerax_path = 'C:/Program Files/ChimeraX 1.8/bin/chimerax.exe'
    ligand_center_output_txt = 'C://Users//klair//Desktop//AutomatedDocking//Automated-Docking//docking//ligand_center.txt'
    ligand_identifier = 'BTN'
    config_file = 'C://Users//klair//Desktop//AutomatedDocking//Automated-Docking//docking//conf.txt'
    output_file = 'C://Users//klair//Desktop//AutomatedDocking//Automated-Docking//docking//docking_output.pdbqt'
    
    # Step 1: Find ligand coordinates
    find_ligand_coordinates(chimerax_path, pdb_file, ligand_identifier, ligand_center_output_txt)
    
    # Step 2: Read coordinates and update config file
    with open(ligand_center_output_txt, 'r') as file:
        lines = file.readlines()
        center_x = lines[0].strip().split('= ')[1]
        center_y = lines[1].strip().split('= ')[1]
        center_z = lines[2].strip().split('= ')[1]

    # Update config file with new coordinates
    with open(config_file, 'w') as file:
        file.write(f"receptor = C://Users//klair//Desktop//AutomatedDocking//Automated-Docking//docking//1stp_DOCK.pdbqt\n")
        file.write(f"ligand = C://Users//klair//Desktop//AutomatedDocking//Automated-Docking//docking//ligand.pdbqt\n")
        file.write(f"center_x = {center_x}\n")
        file.write(f"center_y = {center_y}\n")
        file.write(f"center_z = {center_z}\n")
        file.write(f"size_x = 15\n")
        file.write(f"size_y = 15\n")
        file.write(f"size_z = 15\n")
        file.write(f"cpu = 2\n")
        file.write(f"exhaustiveness = 64\n")
    
    # Step 3: Run Vina
    run_vina(config_file, output_file)
