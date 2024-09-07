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

import os
from django.core.files.storage import FileSystemStorage
from django.conf import settings
from django.core.files.base import ContentFile
from django.shortcuts import render
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import py3Dmol
import subprocess
from .forms import DrugDiscoveryForm
def parse_cohere_output(cohere_output):
    # Define regex patterns to extract required details
    patterns = {
        'nom_ligand': r'Nom du ligand\s*:\s*(.*)',
        'formule_smiles': r'Formule SMILES\s*:\s*(.*)',
        'formule_chimique': r'Formule chimique\s*:\s*(.*)',
        'poids_moleculaire': r'Poids moléculaire\s*:\s*(.*)',
        'configuration': r'Configuration\s*:\s*(.*)',
        'affinite_liaison': r'Affinité de liaison\s*:\s*(.*)',
        'type_liaison': r'Type de liaison\s*:\s*(.*)',
        'site_liaison': r'Site de liaison\s*:\s*(.*)',
    }
    
    # Extract details using regex
    details = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, cohere_output, re.IGNORECASE)
        details[key] = match.group(1).strip() if match else ''

    return details


import numpy as np
from django.http import JsonResponse
from django.shortcuts import render
from django.core.files.base import ContentFile
from django.core.files.storage import FileSystemStorage
import subprocess
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import json
import cohere
from groq import Groq

def drug_discovery_view(request):
    if request.method == 'POST':
        amino_acid_sequence = request.POST.get('amino_acid_sequence')
        virus_name = request.POST.get('virus_name')
        additional_features = request.POST.get('additional_features')

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

        ### Exemple de la maladie d'Alzheimer :
        Examples for answering:
        Question: Pouvez-vous affiner cette structure et ajouter une description détaillée de ses caractéristiques chimiques tel que nom ligand, formule chimique, poids moléculaire,configuration, affinité de liaison, type de liaison interaction avec la protéine cible ainsi que l'importance et la justification de ce choix dans la découverte de médicaments ?
        Answer:
        Nom du ligand: Donepezil
        Structure moléculaire:
                  O=C-N(CH3)-C(CH3)2-C-CN
                  |  |  |   ||   |
                  N-CH3  O  N  CH  N
                  |      |   | \  |
                  C-N-CH3  C  C-C-C-C
                  |    |   |   |   |
                  CH3   H   H   H   H

        Description détaillée:
        Groupe amine tertiaire: Contient un atome d'azote lié à trois groupes alkyle (méthyle et éthyle).
        Ether phénylique: Un atome d'oxygène lié à un groupe phényle (un anneau benzénique substitué par un groupe méthoxy).
        Carbamate: Un groupe ester formé par une liaison entre un atome de carbone carbonyle et un groupe amine.
        Anneau benzénique: Un anneau aromatique à six membres contenant des doubles liaisons conjuguées.
        Groupe méthoxy: Un groupe méthyle lié à un atome d'oxygène.
        Groupe méthyle: Un atome de carbone lié à trois atomes d'hydrogène.
        Détails de la Structure
        Nom du Ligand: Ligand-123
        Formule SMILES: CCCCCCCCCCC(C(=O)O)N
        Formule Chimique: C20H22O4
        Poids Moléculaire: 310.39 g/mol
        Configuration: Configuration tridimensionnelle optimisée avec une géométrie stable.
        Affinité de Liaison
        Affinité de Liaison: 8.5 nM
        Type de Liaison: Liaison covalente
        Interaction avec le Protéine Cible
        Site de Liaison: Site actif de la protéine cible

        Importance et justification de ce choix dans la découverte de médicaments :
        Le donepezil est crucial dans la découverte de médicaments en raison de son rôle en tant qu'inhibiteur de la cholinestérase,
        ce qui augmente les niveaux de neurotransmetteurs dans le cerveau. Utilisé principalement pour traiter la maladie d'Alzheimer,
        il a non seulement démontré une amélioration des fonctions cognitives des patients, mais il a également établi une base pour
        le développement de nouveaux traitements ciblant les voies neurochimiques impliquées dans les maladies neurodégénératives.
        Sa capacité à améliorer les symptômes de la démence a fait de lui un modèle pour la conception de futurs médicaments.
        rappelez-vous, ce n'était qu'un exemple !!.
        ### Demande:
        Pouvez-vous affiner cette structure et ajouter une description détaillée de ses caractéristiques chimiques en utilisant le format suivant :
        - Nom du ligand: 
        - Formule SMILES: 
        - Formule chimique: 
        - Poids moléculaire: 
        - Configuration: 
        - Affinité de liaison: 
        - Type de liaison: 
        - Site de liaison: 
        """

        response = co.generate(
            model="command-r-plus",
            prompt=prompt,
            max_tokens=3000,
            temperature=0.7,
            k=0,
            p=0.75,
            frequency_penalty=0,
            presence_penalty=0,
            stop_sequences=[]
        )

        smiles = smile.strip()
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            return JsonResponse({'error': 'Invalid SMILES string generated. Please try again.'}, status=400)

        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)

        # Save the PDB file in the media directory
        fs = FileSystemStorage(location=settings.MEDIA_ROOT)
        pdb_path = fs.save('ligand_3d.pdb', ContentFile(Chem.MolToPDBBlock(mol)))
        pdb_url = fs.url(pdb_path)

        # Convert the PDB file to PDBQT using Open Babel
        pdbqt_path = os.path.join(settings.MEDIA_ROOT, 'ligand_3d.pdbqt')
        subprocess.run(['obabel', os.path.join(settings.MEDIA_ROOT, pdb_path), '-O', pdbqt_path])

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
        cohere_output = response.generations[0].text.strip()
        parsed_output = parse_cohere_output(cohere_output)

        # Convert NumPy arrays to lists if any
        if isinstance(cent, np.ndarray):
            cent = cent.tolist()

        response_data = {
            'parsed_output': parsed_output,
            'cohere_output': cohere_output,
            'html': html,  # Add additional HTML content if needed
            'pdb_url': '/path/to/pdb',
            'pdbqt_url': '/path/to/pdbqt',
            'cent': 'Some additional data if needed'
        }
        
        return JsonResponse(response_data)
    else:
        return JsonResponse({'error': 'Invalid request method'}, status=400)




import os
import uuid
import requests
import subprocess
from django.conf import settings
from django.shortcuts import render
import os
import uuid
import requests
from django.conf import settings
from django.shortcuts import render

from django.http import JsonResponse
import requests
import uuid
import os
from django.conf import settings
from django.shortcuts import render


def predict_structure(request):
    if request.method == 'POST':
        sequence = request.POST.get('aminoAcidSequence')
        if not sequence:
            return JsonResponse({'error': 'The sequence is required.'})

        try:
            response = requests.post(
                'https://api.esmatlas.com/foldSequence/v1/pdb/',
                headers={'Content-Type': 'text/plain'},
                data=sequence,
                verify=False  # Ignore SSL verification
            )

            protein_id = str(uuid.uuid4())
            pdb_path = os.path.join(settings.MEDIA_ROOT, f'{protein_id}.pdb')

            # Save the PDB file
            with open(pdb_path, 'w') as f:
                f.write(response.text)

            # Load the PDB content
            with open(pdb_path, 'r') as pdb_file:
                pdb_content = pdb_file.read()

            return JsonResponse({'pdb_content': pdb_content})
        except requests.RequestException as e:
            print(f'Request error: {e}')  # Log the exception
            return JsonResponse({'error': 'Error predicting the structure.'})
        except Exception as e:
            print(f'Conversion error: {e}')  # Log the exception
            return JsonResponse({'error': f'Error during conversion: {str(e)}'})

    return render(request, 'finalTemplate2.html')



def convert_pdb_to_mid(pdb_path,mid_path):
    obabel_executable = r'C:\Users\klair\Desktop\AutomatedDocking\Automated-Docking\OpenBabel-2.4.1\babel.exe'
    command = [obabel_executable, pdb_path, '-O', mid_path,'-xr -p 7.4']
    subprocess.run(command, check=True)


    
def convert_mid_to_pdbqt(mid_path,pdbqt_path):
    obabel_executable = r'C:\Users\klair\Desktop\AutomatedDocking\Automated-Docking\OpenBabel-2.4.1\babel.exe'
    command = [obabel_executable, mid_path, '-O', pdbqt_path,'-xr --partialcharge eem']
    subprocess.run(command, check=True)