from django.shortcuts import render


def acc(request):
    return render(request, 'index.html')
