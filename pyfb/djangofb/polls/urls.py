from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('getPositions/',views.getPositions,name='getPositions'),
    path('upload_file/',views.getPositions,name='upload_file')
]
