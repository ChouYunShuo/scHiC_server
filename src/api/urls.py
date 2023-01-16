from django.urls import path
from .views import scHicQueryView, DatasetListAPIView
from . import views
urlpatterns = [
    #path('', main),
    path('query', scHicQueryView.as_view(), name="HiC Contact Map"),
    path('datasets', DatasetListAPIView.as_view()),
]
