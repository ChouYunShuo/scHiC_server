from django.urls import path
from .views import scHicQueryView, DatasetListAPIView, scHicTestView, dataset_retreive_view
from . import views
urlpatterns = [
    #path('', main),
    path('query', scHicQueryView.as_view(), name="HiC Contact Map"),
    #path('datasets', DatasetListAPIView.as_view()),
    path('datasets', dataset_retreive_view),
    path('datasets/<int:pk>/', dataset_retreive_view),
    path('test', scHicTestView.as_view(), name="HiC Contact Map"),
]
