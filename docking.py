import subprocess

def run_vina(config_file, output_file):
    command = [
        'C:/Users/klair/Desktop/DockingJob/Vina/vina.exe',  # Full path to vina executable
        '--config', config_file,
        '--out', output_file
    ]
    subprocess.run(command, check=True)

if __name__ == "__main__":
    config_file = 'C:/Users/klair/Desktop/DockingJob/conf.txt' 
    output_file = 'C:/Users/klair/Desktop/DockingJob/docking_output.pdbqt' 
    run_vina(config_file, output_file)