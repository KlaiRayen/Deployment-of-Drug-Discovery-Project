from django.contrib import admin
from django.urls import path
from docking_interface import views
from django.conf import settings
from django.conf.urls.static import static
from . import views as v
urlpatterns = [
    path('admin/', admin.site.urls),    
    path('upload/', views.upload_files, name='upload'),
    path('upload_pdb/', views.upload_and_visualize_pdb, name='upload_pdb'),
    path('aa', views.upload_and_visualize_pdb, name='upload_and_visualize_pdb'),
    path('', v.acc, name='acc'),
    path('drugDiscovery/', views.drug_discovery_view, name='index'),
    path('predictProtein/', views.predict_structure, name='predict_structure'),
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)