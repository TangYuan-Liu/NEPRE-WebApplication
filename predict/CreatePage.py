import sys
import os
import shutil
def createWaitPage(ID):
    file_path = "./predict/templates/temp/" + str(ID) + "waitpage.html"
    f = open(file_path,'wb')
    savedStdout = sys.stdout
    sys.stdout = f
    
    print'<!DOCTYPE html>'
    print'<html lang="en">'
    print'<head>'
    print'<meta charset="utf-8">'
    print'<title>Nepre-PSPResults</title>'
    print'{% load staticfiles %}'
    print'<meta name="viewport" content="width=device-width, initial-scale=1">'
    print'<meta http-equiv="X-UA-Compatible" content="IE=edge" />'
    print'<link rel="stylesheet" href="{% static "./4/flatly/bootstrap.css" %}" media="screen">'
    print'<link rel="stylesheet" href="{% static "./_assets/css/custom.min.css" %}">'
    print'<script>'
    print'var _gaq = _gaq || [];'
    print"_gaq.push(['_setAccount', 'UA-23019901-1']);"
    print"_gaq.push(['_setDomainName', 'bootswatch.com']);"
    print"_gaq.push(['_setAllowLinker', true]);"
    print"_gaq.push(['_trackPageview']);"
    print'(function() {'
    print"var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;"
    print"ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';"
    print"var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);"
    print'})();'
    print'</script>'
    print'</head>'
    print'<body>'
    print'<div class="navbar navbar-expand-lg fixed-top navbar-dark bg-primary">'
    print'<div class="container">'
    print'<a href="/index/" class="navbar-brand">Nepre</a>'
    print'<button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarResponsive" aria-controls="navbarResponsive" aria-expanded="false" aria-label="Toggle navigation">'
    print'<span class="navbar-toggler-icon"></span>'
    print'</button>'
    print'<div class="collapse navbar-collapse" id="navbarResponsive">'
    print'<ul class="navbar-nav">'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/introduction/">Introduction</a>'
    print'</li>'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/method/">Method</a>'
    print'</li>'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/psp/">Protein Structure Predict</a>'
    print'</li>'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/AminoAcidDistribution/">Amino Acid Distribution</a>'
    print'</li>'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/checkResults/">History</a>'
    print'</li>'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/download/">Download</a>'
    print'</li>'
    print'</ul>'
    print'</div>'
    print'</div>'
    print'</div>'
    print'<div class="page-header" id="banner" style = "margin-left: 80px">'
    print'<div class="row">'
    print'<div class="col-lg-8 col-md-7 col-sm-6">'
    print'<h1>Job Submitted Successfully</h1>'
    print'</div>'
    print'</div>'
    print'</div>'
    print'<div style = "margin-left:80px; margin-top:40px; margin-right:80px">'
    print'<h4 style = "color:cornflowerblue">Processing Scheule</h4>'
    print'<hr></hr>'
    print'</div>'
    print'<div style = "margin-left:120px; margin-top:40px; margin-right:120px">'
    print'<p>Process ID: &nbsp<b>{{ id }}</b></p>'
    print'<p>File Name: &nbsp<b>{{ name }}</b></p>'
    print'<p>'
    print'</div>'
    print'<div style = "margin-left:120px; margin-top:40px; margin-right:80px;text-align:justify">'
    print'<p><b>Please Remember the Job ID carefully, it will be used to get the results.</b></p>' 
    print'<p>Calculation will take some time to finish, please wait.</p>'
    print'<p>Speed: &nbsp<b>150 PDBs/min</b></p>'
    print'</div>'
    
    print'<div style = "margin-left:80px; margin-top:60px; margin-right:80px">'
    print'<div class="alert alert-dismissible alert-info">'
    print'<button type="button" class="close" data-dismiss="alert">&times;</button>'
    print'Click to get<a href=' + '"/checkResults/"' + 'class="alert-link"> results.</a>'
    print'</div>'
    print'</div>'
    print'<script src="{% static "./_vendor/jquery/dist/jquery.min.js" %}"></script>'
    print'<script src="{% static "./_vendor/popper.js/dist/umd/popper.min.js" %}"></script>'
    print'<script src="{% static "./_vendor/bootstrap/dist/js/bootstrap.min.js" %}"></script>'
    print'<script src="{% static "./_assets/js/custom.js" %}"></script>'
    print'<div class="text" style=" text-align:center;margin-top: 120px">Building 9, East Zone, ZPark II, No.10 East Xibeiwang Road, Haidian District, Beijing 100193, China</div>'
    print'<div class="text" style=" text-align:center;">Department of Complex System</div>'
    print'<div class="text" style=" text-align:center;">Email: nepre2018@163.com</div>'
    print'</body>'
    print'</html>'
    
    #reset the stdout
    sys.stdout = savedStdout


"""
def createResultPage(ID, energylist, namelist, folder_path):
    file_path = "./predict/templates/results/" + str(ID) + "results.html"
    decoy_num = len(energylist)
    f = open(file_path,'wb')
    savedStdout = sys.stdout
    sys.stdout = f
    print'<!DOCTYPE html>'
    print'<html lang="en">'
    print'<head>'
    print'<meta charset="utf-8">'
    print'<title>Nepre-PSPResults</title>'
    print'{% load staticfiles %}'
    print'<meta name="viewport" content="width=device-width, initial-scale=1">'
    print'<meta http-equiv="X-UA-Compatible" content="IE=edge" />'
    print'<link rel="stylesheet" href="{% static "./4/flatly/bootstrap.css" %}" media="screen">'
    print'<link rel="stylesheet" href="{% static "./_assets/css/custom.min.css" %}">'
    print'<script>'
    print'var _gaq = _gaq || [];'
    print"_gaq.push(['_setAccount', 'UA-23019901-1']);"
    print"_gaq.push(['_setDomainName', 'bootswatch.com']);"
    print"_gaq.push(['_setAllowLinker', true]);"
    print"_gaq.push(['_trackPageview']);"
    print'(function() {'
    print"var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;"
    print"ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';"
    print"var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);"
    print'})();'
    print'</script>'
    print'</head>'
    print'<body>'
    print'<div class="navbar navbar-expand-lg fixed-top navbar-dark bg-primary">'
    print'<div class="container">'
    print'<a href="/index/" class="navbar-brand">Nepre</a>'
    print'<button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarResponsive" aria-controls="navbarResponsive" aria-expanded="false" aria-label="Toggle navigation">'
    print'<span class="navbar-toggler-icon"></span>'
    print'</button>'
    print'<div class="collapse navbar-collapse" id="navbarResponsive">'
    print'<ul class="navbar-nav">'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/index/">Home</a>'
    print'</li>'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/introduction/">Introduction</a>'
    print'</li>'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/method/">Method</a>'
    print'</li>'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/psp/">Protein Structure Predict</a>'
    print'</li>'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/AminoAcidDistribution/">Amino Acid Distribution</a>'
    print'</li>'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/download/">Download</a>'
    print'</li>'
    print'</ul>'
    print'</div>'
    print'</div>'
    print'</div>'
    print'<div class="page-header" id="banner" style = "margin-left: 80px">'
    print'<div class="row">'
    print'<div class="col-lg-8 col-md-7 col-sm-6">'
    print'<h1>Results</h1>'
    print'</div>'
    print'</div>'
    print'</div>'
    print'<div style="margin-left:80px; margin-right:80px">'
    print'<hr></hr>'
    print'</div>'


    #infornation
    print'<div style = "margin-left:80px; margin-top:40px; margin-right:80px">'
    print'<h4 style = "color:cornflowerblue">Job Information</h4>'
    print'<hr></hr>'
    print'<p>Process ID: {{ id }}</p>'
    print'<p>Total Structures Amount: ' + str(decoy_num) + '</p>'
    print'<p>Submit Time: {{ submit_time }}</p>'
    print'</div>'

    #decoy energy list
    print'<div style = "margin-left:80px; margin-top:20px; margin-right:80px">'
    print'<h4 style = "color:cornflowerblue">Top10 Decoy Energy List</h4>'
    print'<hr></hr>'
    print'</div>'
    print'<div style = "margin-left:400px; margin-top:20px; margin-right:400px">'
    print'<table class="table table-hover">'
    print'<thead>'
    print'<tr>'
    print'<th scope="col">Structure Name</th>'
    print'<th scope="col">Energy</th>'
    print'</tr>'
    print'</thead>'
    print'<tbody>'
    if(decoy_num > 10):
        decoy_dict = {}
        for i in range(decoy_num):
            decoy_dict[namelist[i]] = energylist[i]
        new_dict = sorted(decoy_dict.items(), key = lambda d:d[1])
        for i in range(10):
            print'<tr>'
            print'<th scope="row">' + new_dict[i][0] + '</th>'
            print'<td>' + str(new_dict[i][1]) + '</td>'
            print'</tr>'
    else:
        decoy_dict = {}
        for i in range(decoy_num):
            decoy_dict[namelist[i]] = energylist[i]
        new_dict = sorted(decoy_dict.items(), key = lambda d:d[1])
        for i in range(decoy_num):
            print'<tr>'
            print'<th scope="row">' + new_dict[i][0] + '</th>'
            print'<td>' + str(new_dict[i][1]) + '</td>'
            print'</tr>'
    print'</tbody>'
    print'</table>'
    print'</div>'
    print'<script src="{% static "./_vendor/jquery/dist/jquery.min.js" %}"></script>'
    print'<script src="{% static "./_vendor/popper.js/dist/umd/popper.min.js" %}"></script>'
    print'<script src="{% static "./_vendor/bootstrap/dist/js/bootstrap.min.js" %}"></script>'
    print'<script src="{% static "./_assets/js/custom.js" %}"></script>'
    print'<div class="text" style=" text-align:center;margin-top: 120px">Building 9, East Zone, ZPark II, No.10 East Xibeiwang Road, Haidian District, Beijing 100193, China</div>'
    print'<div class="text" style=" text-align:center;">Department of Complex System</div>'
    print'<div class="text" style=" text-align:center;">Email: nepre2018@163.com</div>'
    print'</body>'
    print'</html>'

    #reset the stdout
    sys.stdout = savedStdout
        
"""


    

def createResultPage(ID, energylist, namelist, folder_path, method):
    dir_path = "./predict/templates/results/" + str(ID) + '/'
    os.mkdir(dir_path)
    file_path = dir_path + str(ID) + "results.html"
    decoy_num = len(energylist)
    f = open(file_path,'wb')
    savedStdout = sys.stdout
    sys.stdout = f

    # sort by energy
    decoy_dict = {}
    for i in range(decoy_num):
        decoy_dict[namelist[i]] = energylist[i]
    new_dict = sorted(decoy_dict.items(), key = lambda d:d[1])

    #copy the lowest decoy here
    target_path = "./uploadfiles/" + str(ID) + '/' + new_dict[0][0]
    desti_path = "./predict/templates/results/" + str(ID) + '/' + new_dict[0][0] 
    shutil.copyfile(target_path, desti_path)

    print'<!DOCTYPE html>'
    print'<html lang="en">'
    print'<head>'
    print'<meta charset="utf-8">'
    print'<meta name="generator" content="Jmol, see http://www.jmol.org">'
    print'<meta name="author" content="siyuan">'
    print'<meta name="keywords" content="generated_by_Jmol_Webexport">'
    print'<title>Nepre-PSPResults</title>'
    print'{% load staticfiles %}'
    print'<meta name="viewport" content="width=device-width, initial-scale=1">'
    print'<meta http-equiv="X-UA-Compatible" content="IE=edge" />'
    print'<link rel="stylesheet" href="{% static "./4/flatly/bootstrap.css" %}" media="screen">'
    print'<link rel="stylesheet" href="{% static "./_assets/css/custom.min.css" %}">'
    print'<script type="text/javascript" src="{% static "./jsmol/JSmol.min.js" %}"></script>'
    print'<script src="{% static "./jsmol/support.js" %}" type="text/javascript"></script>'
    print'<script type="text/javascript" src="{% static "./jsmol/JSmol.min.js" %}"></script>'
    print'<script type="text/javascript">'
    print'var jmolApplet;'
    print'var jmolInfo;'
    print';(function() {'
    print''
    print'var s = document.location.search;'
    print'Jmol.debugCode = (s.indexOf("debugcode") >= 0);'
    print''
    print'jmol_isReady = function(jmolApplet) {'
    print'}'
    print''
    print'jmolInfo = {'
    print'width: "100%",'
    print'height: "100%",'
    print'debug: false,'
    print'color: "white",'
    print'addSelectionOptions: false,'
    print'use: "HTML5",'
    print'coverImage: "@APPLETNAME0@.png",'
    print'coverScript: "",'
    print'deferApplet: false,'
    print'deferUncover: true,'
    print'jarPath: "../../jsmol/java",'
    print'j2sPath: "../../static/jsmol/j2s",'
    print'serverURL: "../../jsmol/php/jsmol.php",'
    print'makeLiveImg:"../../jsmol/j2s/img/play_make_live.jpg",'
    print'jarFile: "JmolAppletSigned0.jar",'
    print'isSigned: true,'
    print''
    print'disableInitialConsole: true,'
    print'readyFunction: jmol_isReady,'
    print'script: ""'
    print'}'
    print''
    print'if (s.indexOf("USE=") >= 0)'
    print'jmolInfo.use = s.split("USE=")[1].split("&")[0]'
    print'else if (s.indexOf("JAVA") >= 0)'
    print'jmolInfo.use = "JAVA"'
    print'else if (s.indexOf("IMAGE") >= 0)'
    print'jmolInfo.use = "IMAGE"'
    print'else if (s.indexOf("NOWEBGL") >= 0)'
    print'jmolInfo.use = "JAVA IMAGE"'
    print'else if (s.indexOf("WEBGL") >= 0)'
    print'jmolInfo.use = "WEBGL HTML5"'
    print'if (s.indexOf("NOWEBGL") >= 0)'
    print'jmolInfo.use = use.replace(/WEBGL/,"")'
    print'var useSignedApplet = (s.indexOf("SIGNED") >= 0);'
    print'if(useSignedApplet && jmolInfo.use == "HTML5") jmolInfo.use = "JAVA";'
    print''
    print'var protocol = window.location.protocol.toLowerCase();'
    print'if (protocol == "file:") {'
    print'jmolInfo.jarPath = "../../jsmol/java";'
    print'jmolInfo.j2sPath = "../../jsmol/j2s";'
    print'jmolInfo.serverURL = "http://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php";'
    print'jmolInfo.makeLiveImg = "../../jsmol/j2s/img/play_make_live.jpg";'
    print'jmolInfo.jarFile = "JmolAppletSigned.jar";'
    print'jmolInfo.isSigned= "true";'
    print'}'
    print'})();'
    print''
    print'$(document).ready(function(){'
    print'var jmolInfo0=jmolInfo;'
    print'jmolInfo0.coverImage="1BYIA.png";'
    print'jmolInfo0.coverScript="javascript revealPopinWidgets(0);";'
    print'jmolInfo0.script="load ' + str(ID) + '"'
    print'$("#Jmol0").html(Jmol.getAppletHtml("jmolApplet0",jmolInfo0));'
    print''
    print'});'
    print'</script>'
    print'<script>'
    print'var _gaq = _gaq || [];'
    print"_gaq.push(['_setAccount', 'UA-23019901-1']);"
    print"_gaq.push(['_setDomainName', 'bootswatch.com']);"
    print"_gaq.push(['_setAllowLinker', true]);"
    print"_gaq.push(['_trackPageview']);"
    print'(function() {'
    print"var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;"
    print"ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';"
    print"var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);"
    print'})();'
    print'</script>'
    print'</head>'
    print'<body>'
    print'<div class="navbar navbar-expand-lg fixed-top navbar-dark bg-primary">'
    print'<div class="container">'
    print'<a href="/index/" class="navbar-brand">Nepre</a>'
    print'<button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarResponsive" aria-controls="navbarResponsive" aria-expanded="false" aria-label="Toggle navigation">'
    print'<span class="navbar-toggler-icon"></span>'
    print'</button>'
    print'<div class="collapse navbar-collapse" id="navbarResponsive">'
    print'<ul class="navbar-nav">'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/introduction/">Introduction</a>'
    print'</li>'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/method/">Method</a>'
    print'</li>'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/psp/">Protein Structure Predict</a>'
    print'</li>'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/AminoAcidDistribution/">Amino Acid Distribution</a>'
    print'</li>'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/checkResults/">History</a>'
    print'</li>'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/download/">Download</a>'
    print'</li>'
    print'</ul>'
    print'</div>'
    print'</div>'
    print'</div>'
    print'<div class="page-header" id="banner" style = "margin-left: 80px">'
    print'<div class="row">'
    print'<div class="col-lg-8 col-md-7 col-sm-6">'
    print'<h1>Results</h1>'
    print'</div>'
    print'</div>'
    print'</div>'
    print'<div style = "margin-left:80px; margin-top:40px; margin-right:80px">'
    print'<h4 style = "color:cornflowerblue">Job Information</h4>'
    print'<hr></hr>'
    print'</div>'
    print'<div style = "margin-left:120px; margin-top:20px; margin-right:120px">'
    print'<p>Process ID:<b>{{ id }}</b></p>'
    print'<p>Method:' + '<b>' + method + '</b>' + '</p>'
    print'<p>Total Structures Amount:' + '<b>' + str(decoy_num) + '</b>'+ '</p>'
    print'<p>Submit Time: <b>{{ submit_time }}</b></p>'
    print'</div>'
    print'<div style = "margin-left:80px; margin-top:40px; margin-right:80px">'
    print'<h4 style = "color:cornflowerblue">Top10 Decoy Energy List</h4>'
    print'<hr></hr>'
    print'</div>'
    print'<div style = "margin-left:400px; margin-top:20px; margin-right:400px">'
    print'<table class="table table-hover">'
    print'<thead>'
    print'<tr>'
    print'<th scope="col">Number</th>'
    print'<th scope="col">Structure Name</th>'
    print'<th scope="col">Energy</th>'
    print'</tr>'
    print'</thead>'
    print'<tbody>'
    if(decoy_num > 10):
        for i in range(10):
            print'<tr>'
            print'<th scope="row">' + str(i) + '</th>'
            print'<td>' + new_dict[i][0] + '</td>'
            print'<td>' + str(new_dict[i][1]) + '</td>'
            print'</tr>'
    else:
        for i in range(decoy_num):
            print'<tr>'
            print'<th scope="row">' + str(i) + '</th>'
            print'<td>' + new_dict[i][0] + '</td>'
            print'<td>' + str(new_dict[i][1]) + '</td>'
            print'</tr>'
    print'</tbody>'
    print'</table>'
    print'</div>'
    print'<div style = "margin-left:120px; margin-top:20px; margin-right:120px">'
    print'<p>Download the <a href="{% static "./download/potential/' + str(ID) + '"%}"> results </a> \
             in text file format.<p>'
    print'</div>'
    print'<div style = "margin-left:80px; margin-top:20px; margin-right:80px">'
    print'<h4 style = "color:cornflowerblue">Highest Probability Structure</h4>'
    print'<hr></hr>'
    print'</div>'
    print'<div style = "margin-left:120px; margin-top:20px; margin-right:120px">'
    print"<p>The prediction is:" + '<b>' + new_dict[0][0] + '</b>' + "</p>"
    print'<div id="section_1BYIA" class="section" style="text-align: center;">'
    print'<table id="scene_1BYIA" class="sceneWrapper floatleft" style="width:300px;margin:auto;">'
    print'<tr>'
    print'<td id="leftJmolCntl0" style="visibility:hidden"></td>'
    print'<td style="width:300px;">'
    print'<div id="Jmol0" class="JmolDiv" style="width:300px; height:300px; background-image:url(1BYIA.png)">'
    print'Jmol0 will appear here.'
    print'</div>'
    print'</td>'
    print'<td id="rightJmolCntl0" style="visibility:hidden"></td>'
    print'</tr>'
    print'<tr>'
    print'</tr>'
    print'</table>'
    print'</div>'
    print'<script src="{% static "./_vendor/jquery/dist/jquery.min.js" %}"></script>'
    print'<script src="{% static "./_vendor/popper.js/dist/umd/popper.min.js" %}"></script>'
    print'<script src="{% static "./_vendor/bootstrap/dist/js/bootstrap.min.js" %}"></script>'
    print'<script src="{% static "./_assets/js/custom.js" %}"></script>'
    print'<div class="text" style=" text-align:center;margin-top: 120px">Building 9, East Zone, ZPark II, No.10 East Xibeiwang Road, Haidian District, Beijing 100193, China</div>'
    print'<div class="text" style=" text-align:center;">Department of Complex System</div>'
    print'<div class="text" style=" text-align:center;">Email: nepre2018@163.com</div>'
    print'</body>'
    print'</html>'

    #reset the stdout
    sys.stdout = savedStdout


def createAAResultPage(amino1,amino2,layer):
    print("&&&&")
    file_path = "./predict/templates/AADistribute/" + amino1 + '-' + amino2 + '-' + layer + ".html"
    f = open(file_path,"wb")
    img_name = amino1 + '-' + amino2 + '-' + layer + ".png"
    savedStdout = sys.stdout
    sys.stdout = f


    print'<!DOCTYPE html>'
    print'<html lang="en">'
    print'<head>'
    print'<meta charset="utf-8">'
    print'<title>Nepre-PSPResults</title>'
    print'{% load staticfiles %}'
    print'<meta name="viewport" content="width=device-width, initial-scale=1">'
    print'<meta http-equiv="X-UA-Compatible" content="IE=edge" />'
    print'<link rel="stylesheet" href="{% static "./4/flatly/bootstrap.css" %}" media="screen">'
    print'<link rel="stylesheet" href="{% static "./_assets/css/custom.min.css" %}">'
    print'<script>'
    print'var _gaq = _gaq || [];'
    print"_gaq.push(['_setAccount', 'UA-23019901-1']);"
    print"_gaq.push(['_setDomainName', 'bootswatch.com']);"
    print"_gaq.push(['_setAllowLinker', true]);"
    print"_gaq.push(['_trackPageview']);"
    print'(function() {'
    print"var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;"
    print"ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';"
    print"var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);"
    print'})();'
    print'</script>'
    print'</head>'
    print'<body>'
    print'<div class="navbar navbar-expand-lg fixed-top navbar-dark bg-primary">'
    print'<div class="container">'
    print'<a href="/index/" class="navbar-brand">Nepre</a>'
    print'<button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarResponsive" aria-controls="navbarResponsive" aria-expanded="false" aria-label="Toggle navigation">'
    print'<span class="navbar-toggler-icon"></span>'
    print'</button>'
    print'<div class="collapse navbar-collapse" id="navbarResponsive">'
    print'<ul class="navbar-nav">'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/index/">Home</a>'
    print'</li>'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/introduction/">Introduction</a>'
    print'</li>'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/method/">Method</a>'
    print'</li>'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/psp/">Protein Structure Predict</a>'
    print'</li>'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/AminoAcidDistribution/">Amino Acid Distribution</a>'
    print'</li>'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/checkResults/">History</a>'
    print'</li>'
    print'<li class="nav-item">'
    print'<a class="nav-link" href="/download/">Download</a>'
    print'</li>'
    print'</ul>'
    print'</div>'
    print'</div>'
    print'</div>'
    print'<div class="page-header" id="banner" style = "margin-left: 80px">'
    print'<div class="row">'
    print'<div class="col-lg-8 col-md-7 col-sm-6">'
    print'<h1>Results</h1>'
    print'</div>'
    print'</div>'
    print'</div>'
    print'<div style="margin-left:80px; margin-right:80px">'
    print'<hr></hr>'
    print'</div>'
    print'<div style="margin-left:200px;margin-right:200px">'
    print'<img src="{% static "./pics/'+ img_name + '"%}"></img>'
    print'</div>'
    #print'<img src="' + img_path + '" alt="pics" style="margin-left:500px">'
    print'<script src="{% static "./_vendor/jquery/dist/jquery.min.js" %}"></script>'
    print'<script src="{% static "./_vendor/popper.js/dist/umd/popper.min.js" %}"></script>'
    print'<script src="{% static "./_vendor/bootstrap/dist/js/bootstrap.min.js" %}"></script>'
    print'<script src="{% static "./_assets/js/custom.js" %}"></script>'
    print'<div class="text" style=" text-align:center;margin-top: 120px">Building 9, East Zone, ZPark II, No.10 East Xibeiwang Road, Haidian District, Beijing 100193, China</div>'
    print'<div class="text" style=" text-align:center;">Department of Complex System</div>'
    print'<div class="text" style=" text-align:center;">Email: nepre2018@163.com</div>'
    print'</body>'
    print'</html>'

    #reset the stdout
    sys.stdout = savedStdout


