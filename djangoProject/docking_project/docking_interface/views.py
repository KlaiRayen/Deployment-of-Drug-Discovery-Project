import os
import subprocess
from django.shortcuts import render
from django.core.files.storage import default_storage
from django.conf import settings
from .forms import UploadFileForm
import re
def find_ligand_coordinates(chimerax_path, pdb_file, ligand_identifier, output_txt_file):
    chimerax_command = (
        f'open {pdb_file}; '
        f'measure center /a:{ligand_identifier};'
        f'close all; exit'
    )
    command = [chimerax_path, '--cmd', chimerax_command]
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print(result.stdout)
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
        'C:/Users/klair/Desktop/DockingJob/Vina/vina.exe',
        '--config', config_file,
        '--out', output_file
    ]
    subprocess.run(command, check=True)

def extract_best_model(pdbqt_file, output_pdb_file):
    with open(pdbqt_file, 'r') as file:
        lines = file.readlines()
    
    best_model_lines = []
    inside_model = False
    model_counter = 0

    for line in lines:
        if line.startswith("MODEL"):
            model_counter += 1
            if model_counter == 2:
                break
            inside_model = True
        elif line.startswith("ENDMDL"):
            inside_model = False

        if inside_model or line.startswith("MODEL 1"):
            best_model_lines.append(line)

    with open(output_pdb_file, 'w') as out_file:
        out_file.writelines(best_model_lines)

def upload_files(request):
    if request.method == 'POST':
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            receptor_file = request.FILES['receptor_file']
            ligand_file = request.FILES['ligand_file']

            receptor_path = os.path.join(settings.MEDIA_ROOT, receptor_file.name)
            ligand_path = os.path.join(settings.MEDIA_ROOT, ligand_file.name)

            with default_storage.open(receptor_path, 'wb+') as destination:
                for chunk in receptor_file.chunks():
                    destination.write(chunk)

            with default_storage.open(ligand_path, 'wb+') as destination:
                for chunk in ligand_file.chunks():
                    destination.write(chunk)

            ligand_center_output_txt = os.path.join(settings.MEDIA_ROOT, 'ligand_center.txt')
            chimerax_path = 'C:/Program Files/ChimeraX 1.8/bin/chimerax.exe'
            ligand_identifier = 'BTN'
            config_file = os.path.join(settings.MEDIA_ROOT, 'conf.txt')
            output_file = os.path.join(settings.MEDIA_ROOT, 'docking_output.pdbqt')
            best_model = os.path.join(settings.MEDIA_ROOT, 'best_model.pdb')

            find_ligand_coordinates(chimerax_path, ligand_path, ligand_identifier, ligand_center_output_txt)

            with open(ligand_center_output_txt, 'r') as file:
                lines = file.readlines()
                center_x = lines[0].strip().split('= ')[1]
                center_y = lines[1].strip().split('= ')[1]
                center_z = lines[2].strip().split('= ')[1]

            with open(config_file, 'w') as file:
                file.write(f"receptor = {receptor_path}\n")
                file.write(f"ligand = {ligand_path}\n")
                file.write(f"center_x = {center_x}\n")
                file.write(f"center_y = {center_y}\n")
                file.write(f"center_z = {center_z}\n")
                file.write(f"size_x = 15\n")
                file.write(f"size_y = 15\n")
                file.write(f"size_z = 15\n")
                file.write(f"cpu = 2\n")
                file.write(f"exhaustiveness = 64\n")

            run_vina(config_file, output_file)
            extract_best_model(output_file, best_model)


            try:
                # Save the receptor file
                with default_storage.open(receptor_path, 'wb+') as destination:
                    for chunk in receptor_file.chunks():
                        destination.write(chunk)

                # Save the ligand file
                with default_storage.open(ligand_path, 'wb+') as destination:
                    for chunk in ligand_file.chunks():
                        destination.write(chunk)

                # Paths for ChimeraX script and output image
                chimerax_script_path = os.path.join(settings.MEDIA_ROOT, 'visualize.cxc')
                output_image_path = os.path.join(settings.MEDIA_ROOT, 'output_image.png')

                # Prepare ChimeraX script
                script_content = prepare_chimerax_script(receptor_path, ligand_path, output_image_path)
                with open(chimerax_script_path, 'w') as script_file:
                    script_file.write(script_content)

                # Path to the ChimeraX executable

                # Generate the URL for the output image
                pdb_url = os.path.join(settings.MEDIA_URL, 'output_image.png')
                return render(request, 'visualize_pdb.html', {'pdb_url': pdb_url})

            except Exception as e:
                logger.error(f'Error during PDB upload and visualization: {e}')
                return render(request, 'upload_pdb.html', {'form': form, 'error': str(e)})

    else:
        form = UploadFileForm()

    return render(request, 'upload_pdb.html', {'form': form})


import os
import subprocess
import logging
from django.conf import settings
from django.core.files.storage import default_storage
from django.shortcuts import render

# Set up logging
logger = logging.getLogger(__name__)

from .forms import UploadFileForm

logger = logging.getLogger(__name__)

def prepare_chimerax_script(receptor_path, ligand_path, output_image_path):
    """Generate ChimeraX script for visualization of receptor and ligand PDB files."""
    # Escape backslashes in the file paths
    receptor_path = receptor_path.replace('\\', '\\\\')
    ligand_path = ligand_path.replace('\\', '\\\\')
    output_image_path = output_image_path.replace('\\', '\\\\')
    
    script_content = f"""
open {receptor_path}
open {ligand_path}
color byhet
view
save {output_image_path} format png
"""
    return script_content



def upload_and_visualize_pdb(request):
    if request.method == 'POST':
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            receptor_file = request.FILES['receptor_file']
            ligand_file = request.FILES['ligand_file']

            # Save the uploaded files
            receptor_path = os.path.join(settings.MEDIA_ROOT, receptor_file.name)
            ligand_path = os.path.join(settings.MEDIA_ROOT, ligand_file.name)

            try:
                # Save the receptor file
                with default_storage.open(receptor_path, 'wb+') as destination:
                    for chunk in receptor_file.chunks():
                        destination.write(chunk)

                # Save the ligand file
                with default_storage.open(ligand_path, 'wb+') as destination:
                    for chunk in ligand_file.chunks():
                        destination.write(chunk)

                # Paths for ChimeraX script and output image
                chimerax_script_path = os.path.join(settings.MEDIA_ROOT, 'visualize.cxc')
                output_image_path = os.path.join(settings.MEDIA_ROOT, 'output_image.png')

                # Prepare ChimeraX script
                script_content = prepare_chimerax_script(receptor_path, ligand_path, output_image_path)
                with open(chimerax_script_path, 'w') as script_file:
                    script_file.write(script_content)

                # Path to the ChimeraX executable
                chimerax_executable = r'C:\\Program Files (x86)\\ChimeraX 1.8\\bin\\ChimeraX.exe'

                # Generate the URL for the output image
                pdb_url = os.path.join(settings.MEDIA_URL, 'output_image.png')
                return render(request, 'visualize_pdb.html', {'pdb_url': pdb_url})

            except Exception as e:
                logger.error(f'Error during PDB upload and visualization: {e}')
                return render(request, 'upload_pdb.html', {'form': form, 'error': str(e)})

    else:
        form = UploadFileForm()

    return render(request, 'upload_pdb.html', {'form': form})



from .forms import DrugDiscoveryForm
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import cohere
from groq import Groq
import py3Dmol

os.environ['GROQ_API_KEY'] = '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'

from django.conf import settings
from django.core.files.base import ContentFile
from django.core.files.storage import FileSystemStorage
from django.shortcuts import render
from .forms import DrugDiscoveryForm
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import py3Dmol
import cohere
from groq import Groq

def extract_ligand_identifier_from_pdb(pdb_file):
    with open(pdb_file, 'r') as file:
        pdb_content = file.read()
    
    # Example: Extracting ligand identifier from PDB
    # This assumes the ligand identifier is part of the HEADER or a specific record type
    match = re.search(r"COMPND\s+1\s+MOLECULE: (.+)", pdb_content)
    if match:
        return match.group(1).strip()
    return "Unknown"

from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def get_center_of_mass(mol):
    # Add hydrogens if necessary
    mol = Chem.AddHs(mol)
    # Compute 3D coordinates if not present
    if not mol.GetNumConformers():
        AllChem.EmbedMolecule(mol)
    conf = mol.GetConformer()
    # Get atomic coordinates
    pos = conf.GetPositions()
    # Compute center of mass
    center_of_mass = np.mean(pos, axis=0)
    return center_of_mass

def drug_discovery_view(request):
    if request.method == 'POST':
        form = DrugDiscoveryForm(request.POST)
        if form.is_valid():
            amino_acid_sequence = form.cleaned_data['amino_acid_sequence']
            virus_name = form.cleaned_data['virus_name']
            additional_features = form.cleaned_data['additional_features']

            llama3_user_prompt = f"""
                You are a researcher working on chemical structures. Provide only and without any description a SMILES structure for a ligand that does not already exist to target a virus for use in drug discovery. The ligand should only contain the following features:
                - benzene rings
                - carboxylic acid groups
                - amine groups
                - side chains
                Give only and only the SMILES representation. No descriptions or details or introduction to what are you going to give.
                The virus details are as follows:
                Name: {virus_name}
                Amino Acid Sequence:
                {amino_acid_sequence}
                Additional Characteristics:
                {additional_features}
            """

            client = Groq()
            completion = client.chat.completions.create(
                model="llama3-70b-8192",
                messages=[{"role": "user", "content": llama3_user_prompt}],
                temperature=1,
                max_tokens=1024,
                top_p=1,
                stream=True,
                stop=None,
            )

            smile = ""
            for chunk in completion:
                smile += chunk.choices[0].delta.content or ""

            # Cohere part
            co = cohere.Client(api_key="5D2ZAmR16ezQCY893TAP8txMO24yEFMbK8wBbAyA")
            prompt = f"""
            Vous êtes un chercheur spécialisé dans la découverte de médicaments. Voici une structure SMILE générée par un modèle :
            {smile}
            Pouvez-vous affiner cette structure et ajouter des détails supplémentaires, notamment des propriétés pharmacocinétiques potentielles et une description détaillée ?
            """

            response = co.generate(
                model="command-r-plus",
                prompt=prompt,
                max_tokens=300,
                temperature=0.7,
                k=0,
                p=0.75,
                frequency_penalty=0,
                presence_penalty=0,
                stop_sequences=[]
            )

            def validate_smiles(smiles):
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is None:
                        return False
                    return True
                except:
                    return False

            smiles = smile.strip()
            if not validate_smiles(smiles):
                return render(request, 'drug_discovery/invalid_smiles.html')

            mol = Chem.MolFromSmiles(smiles)
            img = Draw.MolToImage(mol)

            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)

            # Save the PDB file in the media directory
            fs = FileSystemStorage(location=settings.MEDIA_ROOT)
            pdb_path = fs.save('ligand_3d.pdb', ContentFile(Chem.MolToPDBBlock(mol)))
            pdb_url = fs.url(pdb_path)

            # Extract ligand identifier from PDB
            ligand_identifier = extract_ligand_identifier_from_pdb(os.path.join(settings.MEDIA_ROOT, 'ligand_3d.pdb'))

            # Prepare 3D visualization
            view = py3Dmol.view(width=800, height=400)
            with open(settings.MEDIA_ROOT + '/ligand_3d.pdb', 'r') as f:
                pdb_data = f.read()
            view.addModel(pdb_data, 'pdb')
            view.setStyle({'stick': {}})
            view.setBackgroundColor('white')
            view.zoomTo()
            html = view._make_html()

            cent = get_center_of_mass(mol)

            context = {
                'form': form,
                'smiles': smile,
                'cohere_output': response.generations[0].text,
                'img': img,
                'html': html,
                'pdb_url': pdb_url,  # Add the URL for downloading the PDB file
                'ligand_identifier': ligand_identifier,
                'cent': cent,
            }
            return render(request, 'drug_discovery/result.html', context)
    else:
        form = DrugDiscoveryForm()
    return render(request, 'drug_discovery/index.html', {'form': form})





import os
import uuid
import requests
import subprocess
from django.conf import settings
from django.shortcuts import render

def predict_structure(request):
    if request.method == 'POST':
        sequence = request.POST.get('sequence')
        if not sequence:
            return render(request, 'proteins/predict.html', {'error': 'La séquence est requise.'})

        try:
            response = requests.post(
                'https://api.esmatlas.com/foldSequence/v1/pdb/',
                headers={'Content-Type': 'text/plain'},
                data=sequence,
                verify=False  # Ignore SSL verification
            )

            protein_id = str(uuid.uuid4())
            pdb_path = os.path.join(settings.MEDIA_ROOT, f'{protein_id}.pdb')
            pdbqt_path = os.path.join(settings.MEDIA_ROOT, f'{protein_id}.pdbqt')

            # Save the PDB file
            with open(pdb_path, 'w') as f:
                f.write(response.text)

            # Convert PDB to PDBQT using Open Babel command line
            convert_pdb_to_pdbqt(pdb_path, pdbqt_path)

            return render(request, 'proteins/predict.html', {'message': 'La prédiction de la structure a été sauvegardée.', 'file': f'{protein_id}.pdbqt', 'id': protein_id})
        except requests.RequestException as e:
            return render(request, 'proteins/predict.html', {'error': 'Erreur lors de la prédiction de la structure.'})
        except Exception as e:
            return render(request, 'proteins/predict.html', {'error': f'Erreur lors de la conversion: {str(e)}'})

    return render(request, 'proteins/predict.html')

def convert_pdb_to_pdbqt(pdb_path, pdbqt_path):
    # Path to Open Babel executable
    obabel_executable = r'C:\Users\klair\Desktop\AutomatedDocking\Automated-Docking\OpenBabel-2.4.1\babel.exe'
    command = [obabel_executable, pdb_path, '-O', pdbqt_path]
    subprocess.run(command, check=True)

