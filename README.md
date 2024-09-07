# Deployment of Drug Discovery Project

## Overview

In this project, you'll explore the full pipeline of a drug discovery process, focusing on the interaction between viral proteins and potential therapeutic compounds. The key steps include:

1. **3D Protein Structure Modeling**  
   Using ESMFold, we accurately predict the 3D structures of viral proteins. This allows us to identify crucial binding sites essential for the drug development process. 🧬🔍

2. **Drug Ligand Generation**  
   Our multimodal large language model fine-tunes and optimizes drug ligands to precisely target the identified viral proteins.💡💊

3. **Binding Site Identification**  
   Leveraging AutoDock-Vina, we identify specific regions on proteins where ligands can bind, enhancing the precision and effectiveness of potential treatments. 🎯🔗

## Demo  
Check out the demo here: [LinkedIn Demo](https://www.linkedin.com/posts/rayen-klai-910086244_summercamp2024-biology-biotech-activity-7233100823147261952-Xxuu?utm_source=share&utm_medium=member_desktop)  
Feel free to follow me for updates!

## How to Run the Project

To run this Django project locally, follow these steps:

1. Clone the repository to your local machine:
   ```bash
   git clone <repository_url>
   ```

2. Navigate to the project directory:
   ```bash
   cd <project_directory>
   ```

3. Install the required dependencies:
   ```bash
   pip install -r requirements.txt
   ```

4. Apply the migrations:
   ```bash
   python manage.py migrate
   ```

5. Start the Django development server:
   ```bash
   python manage.py runserver
   ```

6. Access the application at:
   ```
   http://127.0.0.1:8000/
   ```

Feel free to reach out if you encounter any issues or have feedback!
