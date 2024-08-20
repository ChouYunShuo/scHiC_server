from django.urls import path
from .views import scHicQueryView, chromLenView, dataset_retreive_view,embeddingView, spatialView, trackView, geneExprView, metaView, session_retreive_view, session_upload_view
from . import views
urlpatterns = [
    #path('', main),
    path('query', scHicQueryView.as_view(), name="HiC Contact Map"),
    path('datasets', dataset_retreive_view),
    path('datasets/<int:pk>/', dataset_retreive_view),
    path('session/<str:uuid>/', session_retreive_view),
    path('chromlens', chromLenView),
    path('embed', embeddingView),
    path('spatial', spatialView),
    path('meta', metaView),
    path('gene_expr', geneExprView),
    path('track', trackView),
    path('session_upload', session_upload_view)
]
