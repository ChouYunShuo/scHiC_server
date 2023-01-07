from django.urls import path
from .views import scHicView

urlpatterns = [
    #path('', main),
    path('query', scHicView.as_view(), name="HiC Contact Map")
]
