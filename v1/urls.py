"""v1 URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.11/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf.urls import url
from django.contrib import admin
from predict import views as predictViews

urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^index/', predictViews.index),
    url(r'^introduction/', predictViews.introduction),
    url(r'^psp/', predictViews.psp),
    url(r'^nepre_r/',predictViews.nepre_r),
    url(r'^nepre_f/',predictViews.nepre_f),
    url(r'^checkResults/',predictViews.checkResults),
    url(r'^getResultsPage/',predictViews.getResultsPage),
    url(r'^AminoAcidDistribution/',predictViews.AADistribute),
    url(r'^getAADistribute/',predictViews.getAADistribute),
    url(r'^download/',predictViews.download),
    url(r'^method/', predictViews.method),
    url(r'^primary/', predictViews.primary),
]
