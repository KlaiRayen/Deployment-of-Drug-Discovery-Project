from django.shortcuts import render

print("Packages installed successfully!")


def acc(request):
    return render(request, 'index.html')

from docking_interface.forms import DrugDiscoveryForm

def interface(request):
    form = DrugDiscoveryForm(request.POST)
    return render(request, 'interface/interface1.html',{'form' : form})

def pipeline(request):
    return render(request, 'pipeline.html')