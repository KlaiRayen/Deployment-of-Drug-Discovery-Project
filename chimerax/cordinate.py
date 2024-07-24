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

if __name__ == "__main__":
    pdb_file = 'C:/Users/klair/Desktop/DockingJob/1stp.pdb' 
    chimerax_path = 'C:/Program Files/ChimeraX 1.8/bin/chimerax.exe'  
    ligand_center_output_txt = 'C:/Users/klair/Desktop/DockingJob/ligand_center.txt'
    ligand_identifier = 'BTN'

    find_ligand_coordinates(chimerax_path, pdb_file, ligand_identifier, ligand_center_output_txt)