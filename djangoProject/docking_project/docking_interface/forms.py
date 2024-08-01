from django import forms

class UploadFileForm(forms.Form):
    receptor_file = forms.FileField(label='Receptor File')
    ligand_file = forms.FileField(label='Ligand File')


class DrugDiscoveryForm(forms.Form):
    amino_acid_sequence = forms.CharField(widget=forms.Textarea, label="Enter the amino acid sequence of the virus")
    virus_name = forms.CharField(widget=forms.Textarea, label="Enter the name of the virus")
    additional_features = forms.CharField(widget=forms.Textarea, label="Enter other characteristics of the virus")