# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.shortcuts import render
import os
import thread
import time
import zipfile
import SingleStructureCalculate as sc 
from CreatePage import createWaitPage, createResultPage, createAAResultPage,generate_3d_view
from django.http import HttpResponse
import MakePlot
import sys
sys.path.append("./predict/Radius")
sys.path.append("./predict/Cutoff")
import Nepre_R
import Nepre_F
import smtplib
from email.mime.text import MIMEText
from email.header import Header
from email.mime.multipart import MIMEMultipart
# Create your views here.

def index(request):
    pass
    return render(request, "index.html")

def introduction(request):
    pass
    return render(request, "introduction.html")

def method(request):
    pass
    return render(request, "method.html")

def psp(request):
    pass
    return render(request, "psp.html")


def send_email(SMTP_host, from_account, from_passwd, to_account, subject, content):

    email_client = smtplib.SMTP(SMTP_host)
    email_client.login(from_account, from_passwd)
    # create msg
    msg = MIMEText(content, 'plain', 'utf-8')
    msg['Subject'] = Header(subject, 'utf-8')
    msg['From'] = from_account
    msg['To'] = to_account
    email_client.sendmail(from_account, to_account, msg.as_string())
    email_client.quit()
    print("Done")

def send_notice(text,to_account):

    host = "smtp.163.com"
    username = "nepre2018@163.com"
    passwd = "nepre2018liulab"
    subject = "NEPRE Science"
    send_email(host, username, passwd, to_account, subject, text)



def calcEnergy_r(request,cur_time,email):

    if request.method == "POST":
        file_obj = request.FILES.get("up_file")
        cur = str(cur_time)
        os.mkdir("./uploadfiles/" + cur)
        f1 = open("./uploadfiles/" + cur + '/' + file_obj.name, "wb")
        base_path = "./uploadfiles/" + cur + '/'

        for i in file_obj.chunks():
            f1.write(i)
        f1.close()
        print(file_obj.name)
        
        # unzip the files
        namelist = []
        zf = zipfile.ZipFile(base_path+file_obj.name)
        for name in zf.namelist():
            namelist.append(name)
            zf.extract(name,path = base_path)
        # start to calculate the energy
        energy_list = []
        Matrix = Nepre_R.load_EnergyMatrix("./predict/Radius/radius.npy")
        for name in namelist:
            decoy_path = base_path + name
            df = open(decoy_path)
            E = Nepre_R.calculate_Energy(df,Matrix,"./predict/Radius/mean_radius.txt")
            energy_list.append(E)
        folder_path = "./results"
        # write the results in a text file
        res_file = open("./predict/static/download/potential/" + cur, "wb")
        res_file.write("Nepre-R Potential Created by LiuLab of Beijing Computation Science Research Center.")
        res_file.write("\n")
        for i in range(len(energy_list)):
            res_file.write(namelist[i] + "  " + str(energy_list[i]) + "\n")
        res_file.close()
        createResultPage(cur_time,energy_list,namelist,folder_path,"Nepre-R")

    print "Cal Finish:",time.time()
    print energy_list
    
    text = "Your job submitted in NEPRE-WebApplication has completed." + '\n' + "ID: " + cur
    send_notice(text,email)


def calcEnergy_f(request,cur_time,email):

    if request.method == "POST":
        file_obj = request.FILES.get("up_file")
        cur = str(cur_time)
        os.mkdir("./uploadfiles/" + cur)
        f1 = open("./uploadfiles/" + cur + '/' + file_obj.name, "wb")
        base_path = "./uploadfiles/" + cur + '/'
        cutoff = int(request.POST["cutoff"])

        for i in file_obj.chunks():
            f1.write(i)
        f1.close()
        print(file_obj.name)
        
        # unzip the files
        namelist = []
        zf = zipfile.ZipFile(base_path+file_obj.name)
        for name in zf.namelist():
            namelist.append(name)
            zf.extract(name,path = base_path)
        # start to calculate the energy
        energy_list = []
        Matrix = Nepre_F.load_EnergyMatrix(cutoff,"./predict/Cutoff/")
        for name in namelist:
            decoy_path = base_path + name
            df = open(decoy_path)
            E = Nepre_F.calculate_Energy(df,Matrix,cutoff)
            energy_list.append(E)
        folder_path = "./results"

        # write the results in a text file
        res_file = open("./predict/static/download/potential/" + cur, "wb")
        res_file.write("Nepre-F Potential Created by LiuLab of Beijing Computation Science Research Center.")
        res_file.write("\n")
        res_file.write("Cutoff: " + str(cutoff) + " angstrom")
        res_file.write("\n")
        for i in range(len(energy_list)):
            res_file.write(namelist[i] + "  " + str(energy_list[i]) + "\n")
        res_file.close()
        createResultPage(cur_time,energy_list,namelist,folder_path,"Nepre-F")

    print "Cal Finish:",time.time()
    print energy_list
    text = "Your job submitted in NEPRE-WebApplication has completed." + '\n' + "ID: " + cur
    send_notice(text,email)

def nepre_r(request):
    name = ""
    if request.method == "POST":
        file_obj = request.FILES.get("up_file")
        name = file_obj.name
        email = request.POST["to_email"]
    cur_time = int(time.time())
    try:
        thread.start_new_thread(calcEnergy_r, (request,cur_time,email,))
    except:
        return render(request,"error.html")
    print "psp finish:", time.ctime()
    createWaitPage(cur_time)
    context = {}
    context["id"] = str(cur_time)
    context["name"] = name
    path = "./temp/" + str(cur_time) + "waitpage.html"
    return render(request, path, context)    

def nepre_f(request):
    name = ""
    if request.method == "POST":
        file_obj = request.FILES.get("up_file")
        name = file_obj.name
    cur_time = int(time.time())
    try:
        thread.start_new_thread(calcEnergy_f, (request,cur_time,))
    except:
        return render(request,"error.html")
    print "psp finish:", time.ctime()
    createWaitPage(cur_time)
    context = {}
    context["id"] = str(cur_time)
    context["name"] = name
    path = "./temp/" + str(cur_time) + "waitpage.html"
    return render(request, path, context)



def checkResults(request):
    pass
    return render(request, "checkResults.html")


def getResultsPage(request):
    
    if request.method == "POST":
        Job_id = request.POST["job_id"]
        print request.POST
    if request.method == "GET":
        url = request.path
        url = url.encode()
        print("This is URL")
        print url
        
        struc = url[-10:]
        #print struc
        
        struc_dir = "./predict/templates/results/" + struc + "/"
        file_list = os.listdir(struc_dir)
        file_name = None
        for name in file_list:
            if(name[-3:] == "pdb"):
                file_name = name
                break
        struc_path = struc_dir + file_name
        print struc_path
        
        #f = open("./predict/templates/results/5mlm.pdb",'rb')
        f = open(struc_path,'rb')
        return HttpResponse(f,content_type="text/pdb")

    #job_id = request.POST["job_id"]
    html_path = "./results/" + Job_id + '/' + Job_id + "results.html"
    context = {}
    context["id"] = Job_id
    context["submit_time"] = time.ctime(int(Job_id))
    return render(request, html_path, context)


def AADistribute(request):
    pass
    return render(request, "AminoAcidDistribution.html") 
    

def getAADistribute(request):
    if request.method == "POST":
        #layer = request.POST["layer"]
        center = request.POST["center"]
        surround = request.POST["surround"]
        #print layer,center,surround
        
        #os.system("python ./predict/MakePlot.py " + center + " " + surround + " " + '0')
        #MakePlot.contourmap(center,surround,0)
        
        #createAAResultPage(center,surround,layer)
        generate_3d_view(center,surround)
        html_path = "./AADistribute/" + center + '-' + surround + ".html"
        
        return render(request, html_path)

    if request.method == "GET":
        url = request.path
        url = url.encode()
        print("This is URL")
        print url
        data = url[-11:]
        
        data_path = "./predict/AA3D_view/" + data
        f = open(data_path,'rb')
        return HttpResponse(f,content_type="text/pdb")


def download(request):
    pass
    return render(request, "download.html")
