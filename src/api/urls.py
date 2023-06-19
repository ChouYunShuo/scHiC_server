from django.urls import path
from .views import scHicQueryView, chromLenView, dataset_retreive_view,embeddingView
from . import views
urlpatterns = [
    #path('', main),
    path('query', scHicQueryView.as_view(), name="HiC Contact Map"),
    path('datasets', dataset_retreive_view),
    path('datasets/<int:pk>/', dataset_retreive_view),
    path('chromlens', chromLenView),
    path('embed', embeddingView)
    #path('test', scHicTestView.as_view(), name="HiC Contact Map"),
]
