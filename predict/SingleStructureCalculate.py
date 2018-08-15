import os
import math
import numpy as np
import AminoAcid as AA
import gc
import sys
import math
import csv
import matplotlib.pyplot as plt
from matplotlib.pyplot import savefig



aminoDict={"ALA":[0]*15,"VAL":[0]*15,"LEU":[0]*15,"ILE":[0]*15,"PHE":[0]*15,\
           "TRP":[0]*15,"MET":[0]*15,"PRO":[0]*15,"GLY":[0]*15,"SER":[0]*15,\
           "THR":[0]*15,"CYS":[0]*15,"TYR":[0]*15,"ASN":[0]*15,"GLN":[0]*15,\
           "HIS":[0]*15,"LYS":[0]*15,"ARG":[0]*15,"ASP":[0]*15,"GLU":[0]*15,}

"""
cdDict={"ALA":{},"VAL":{},"LEU":{},"ILE":{},"PHE":{},\
        "TRP":{},"MET":{},"PRO":{},"GLY":{},"SER":{},\
        "THR":{},"CYS":{},"TYR":{},"ASN":{},"GLN":{},\
        "HIS":{},"LYS":{},"ARG":{},"ASP":{},"GLU":{},}
"""
def LoadRadius(path):
    radiusDict = {"ALA":0,"VAL":0,"LEU":0,"ILE":0,"PHE":0,\
                  "TRP":0,"MET":0,"PRO":0,"GLY":0,"SER":0,\
                  "THR":0,"CYS":0,"TYR":0,"ASN":0,"GLN":0,\
                  "HIS":0,"LYS":0,"ARG":0,"ASP":0,"GLU":0,}

    f = open(path)
    for line in f.readlines():
        temp = line.strip().split()
        if(temp[0] != "Name"):
            radiusDict[temp[0]] = float(temp[3])
    return radiusDict


def pearson(rmsd,energy):
    size = np.shape(rmsd)[0]
    x=np.empty(shape=[2,size])
    for i in range(size):
        x[0][i]=rmsd[i]
    for j in range(size):
        x[1][j]=energy[j]
    y=np.corrcoef(x)
    return y[0][1]


def loadModel(path,cdDict):
    List = cdDict.keys()
    List.sort()
    
    f1 = open(path,"rb")
    for amino1 in List:
        for amino2 in List:
            cdDict[amino1][amino2] = np.load(f1)
    f1.close()
    return cdDict

def StringSolve(s):
    blank = 0
    content = 0
    count = 0
    for i in range(len(s)):
        if(s[i] == ' '):
            blank = i
            break
    DecoyName = s[:blank]
    temp = s[blank:]
    for i in range(len(temp)):
        if(temp[i] != ' ' and count == 0):
            content = i
            count += 1
        else:
            if(temp[i] != ' '):
                count += 1
            else:
                if(temp[i] == ' ' and count != 0):
                    break
    RMSD = temp[content:content+count]
    return DecoyName,RMSD

def ExtractData(line):
    String = line.strip()
    start = None
    end = 0
    flag = 0
    ElementList = []
    AAtype = String[17:20]
    AANUMBER = String[24:27]
    for i in range(len(String)):
        if(i == len(String) - 1):
            if(flag == 1):
                ElementList.append(String[start:i+1])
        if(String[i] != ' ' and flag == 0):
            start = i
            flag = 1
        else:
            if(String[i] == ' ' and flag == 1):
                ElementList.append(String[start:i])
                flag = 0
                start = None
    if(len(ElementList) < 10 or len(ElementList) == 10):
        if(len(ElementList) != 1):
            #This line has problem
            templist = ElementList
            flag = None
            temp = []
            Seperate = []
            for i in range(len(ElementList)):
                if(len(ElementList[i]) > 10):
                    flag = i
                    break
            if(flag == None):
                flag = 2
            if(flag == 2):
                Seperate.append(ElementList[2][0:3])
                Seperate.append(ElementList[2][3:])
                for i in range(len(templist)):
                    if(i == flag):
                        ElementList.extend(Seperate)
                    else:
                        ElementList.append(templist[i])
            else:
                element = ElementList[flag]
                
                
                if(element[0] != '-'):
                    temp.append(0)
                
                for i in range(len(element)):
                    if(element[i] == '-'):
                        temp.append(i)
                
                    
                for i in range(len(temp)):
                    if(i == len(temp) - 1):
                        Seperate.append(element[temp[i]:len(element)])
                    else:
                        Seperate.append(element[temp[i]:temp[i+1]])
                ElementList = []
                for i in range(len(templist)):
                    if(i == flag):
                        ElementList.extend(Seperate)
                    else:
                        ElementList.append(templist[i])
    return ElementList,AAtype,AANUMBER


def SingleStructure(decoyname,DecoyPath,model_path,radius_path):
    
    cdDict={"ALA":{},"VAL":{},"LEU":{},"ILE":{},"PHE":{},\
        "TRP":{},"MET":{},"PRO":{},"GLY":{},"SER":{},\
        "THR":{},"CYS":{},"TYR":{},"ASN":{},"GLN":{},\
        "HIS":{},"LYS":{},"ARG":{},"ASP":{},"GLU":{},}

    cdDict = loadModel(model_path,cdDict)
    radiusDict = LoadRadius(radius_path)
    

    if(decoyname == "native.pdb"):   
        xcor = 6
        ycor = 7
        zcor = 8
        AAnum = 5
    else:
        xcor = 5
        ycor = 6
        zcor = 7
        AAnum = 4
    
   
    df = open(DecoyPath)
    #CurrentAAName = None
    CurrentAANitrogen = None
    CurrentAACA = None
    CurrentAANumber = None
    EachAA = []
    
    CurrentAA = None 


    for line in df.readlines():
        #print line
        Element,AAtype,AANUMBER = ExtractData(line)
        
        if(Element[0] != "ATOM"):
            CurrentAA.CalculateCenter()
            CurrentAA.InputCAN(CurrentAANitrogen,CurrentAACA)
            EachAA.append(CurrentAA)
            continue
        
        if(Element[2] == "H"):
            continue
        
        if(AAtype not in cdDict):
            continue
        if(CurrentAA == None):
            #print("First object establised")
            #CurrentAAName = Element[3]
            CurrentAA = AA.AminoAcid(AAtype)
            CurrentAANumber = AANUMBER
            if(Element[2] == "N" or Element[2] == "CA"):
                if(line[16] == "B"):
                    continue
                if(Element[2] == "N"):
                    CurrentAANitrogen = np.array([float(Element[xcor]),float(Element[ycor]),float(Element[zcor])])
                else:
                    CurrentAACA = np.array([float(Element[xcor]),float(Element[ycor]),float(Element[zcor])])
            if(AAtype == "GLY" or Element[2] not in {"N","CA","C","O","O1","02"}):
                if(line[16] != " "):
                    #If cases like "AASN or BASN" appears, we only add A 
                    if(line[16] == "A" and line[15] == "1"):
                        CurrentAA.SumCenters(float(Element[xcor]),float(Element[ycor]),float(Element[zcor]))
                    else:
                        continue
                else:
                    CurrentAA.SumCenters(float(Element[xcor]),float(Element[ycor]),float(Element[zcor]))
        else:
            #If another amino acid begins
            if(AANUMBER != CurrentAANumber):
                #print CurrentAA.AminoAcidAmount
                #print CurrentAAName,Element[3]
                state = CurrentAA.CalculateCenter()
                if(state == False):
                    CurrentAA = AA.AminoAcid(AAtype)
                    CurrentAANumber = AANUMBER
                    continue
                CurrentAA.InputCAN(CurrentAANitrogen,CurrentAACA)
                #print sys.getrefcount(CurrentAA)
                EachAA.append(CurrentAA)
                del CurrentAA
                CurrentAA = AA.AminoAcid(AAtype)
                #print sys.getrefcount(CurrentAA)
                CurrentAANumber = AANUMBER
                if(Element[2] == "N" or Element[2] == "CA"):
                    if(line[16] == "B"):
                        continue
                    if(Element[2] == "N"):
                        CurrentAANitrogen = np.array([float(Element[xcor]),float(Element[ycor]),float(Element[zcor])])
                    else:
                        CurrentAACA = np.array([float(Element[xcor]),float(Element[ycor]),float(Element[zcor])])
                if(AAtype == "GLY" or Element[2] not in {"N","CA","C","O","O1","02"}):
                    if(line[16] != " "):
                    #If cases like "AASN or BASN" appears, we only add A 
                        if(line[16] == "A" and line[15] == "1"):
                            CurrentAA.SumCenters(float(Element[xcor]),float(Element[ycor]),float(Element[zcor]))
                        else:
                            continue
                    else:
                        CurrentAA.SumCenters(float(Element[xcor]),float(Element[ycor]),float(Element[zcor]))
            #If still the same amino acid
            else:
                if(Element[2] == "N" or Element[2] == "CA"):
                    if(line[16] == "B"):
                        continue
                    if(Element[2] == "N"):
                        CurrentAANitrogen = np.array([float(Element[xcor]),float(Element[ycor]),float(Element[zcor])])
                    else:
                        CurrentAACA = np.array([float(Element[xcor]),float(Element[ycor]),float(Element[zcor])])
                if(AAtype == "GLY" or Element[2] not in {"N","CA","C","O","O1","02"}):
                    if(line[16] != " "):
                    #If cases like "AASN or BASN" appears, we only add A 
                        if(line[16] == "A" and line[15] == "1"):
                            CurrentAA.SumCenters(float(Element[xcor]),float(Element[ycor]),float(Element[zcor]))
                        else:
                            continue
                    else:
                        CurrentAA.SumCenters(float(Element[xcor]),float(Element[ycor]),float(Element[zcor]))
        
        
    del CurrentAA#Free the current object.                
    #Scan over. Each amino acid is stored as an object in EachAA. Next step is to calculate the energy, results will be saved in EnergyList. 
    E = 0 #Store the energy
    Time = 0
    for m in range(len(EachAA)):
        #Establish axis first    
        EachAA[m].EstablishCoordinate()
        for n in range(len(EachAA)):
            if(m == n):
                continue
            else:
                dis = EachAA[m].DistanceBetweenAA(EachAA[n].center)
                radiusSum = radiusDict[EachAA[m].name] + radiusDict[EachAA[n].name]
                if(dis <= radiusSum):#If the distance between two amino acid less than 10, we believe the two amino acid have interaction
                    #print EachAA[m].ChangeCoordinate(EachAA[n].center)   
                    rho,theta,phi = EachAA[m].ChangeCoordinate(EachAA[n].center)
                    theta = min(int(math.floor(theta*20/np.pi)),19)
                    phi = min(int(math.floor(phi*10/np.pi) + 10),19)
                    
                    #print EachAA[m].name,EachAA[n].name
                    E += cdDict[EachAA[m].name][EachAA[n].name][theta][phi] / rho 
                    Time += 1
                    
    return E,Time


if __name__ == "__main__":
    """
    Command Example:
    python SingleStructureCalculate 1BYIA 3DRobot_set P
   
    Paramters:
    P  ----- Show Plot Only
    PS ----- Show Plot and Save Results
    S  ----- Only Save Results(Figure is saved in folder "./PicResults")
    """
    args = sys.argv[1:]
    structurename = args[0]
    options = args[2]
    TestDataset = args[1]
    #loadModel()
    path = "../../" + TestDataset  + '/' + structurename + "/list.txt"
    DecoyList = []
    DecoyRMSD = []
    ResultList = []
    TimeList = []
    f = open(path)
    #Read the nature structure and decoy
    for line in f.readlines():
        if line[:4] != "NAME":
            line = line.strip()
            dname,drmsd = StringSolve(line)
            DecoyList.append(dname)
            DecoyRMSD.append(float(drmsd))
    f.close()
    
    for i in range(len(DecoyList)):
        decoypath = "../../" + TestDataset + '/' + structurename + '/' + DecoyList[i]
        E,Time = SingleStructure(DecoyList[i],decoypath)
        ResultList.append(E)
        #TimeList.append(Time)

    if(options == 'P'):
        plt.figure()
        plt.scatter(DecoyRMSD,ResultList,color = "red")
        plt.xlabel("RMSD")
        plt.ylabel("Energy")
        plt.title(structurename)
        plt.show()
        plt.close()

    if(options == "PS" ):
        plt.figure()
        plt.scatter(DecoyRMSD,ResultList,color = "red")
        plt.xlabel("RMSD")
        plt.ylabel("Energy")
        plt.title(structurename)
        savefig("./Results/" + TestDataset + "/Pics/" + structurename + '.png')
        plt.show(block = False)
        plt.close()
        EFile = open("./Results/" + TestDataset + "/Energy/Energy_" + structurename + ".txt","a")
        for i in range(len(ResultList)):
            EFile.write(str(ResultList[i]))
            EFile.write("\t")
            EFile.write(str(DecoyRMSD[i]))
            EFile.write("\n")
        EFile.close()
        pea = pearson(DecoyRMSD,ResultList)
        ENERGY = np.asarray(ResultList)
        INDEX = np.argsort(ENERGY)
        NATIVEINDEX = None
        for i in range(10):
            if(DecoyList[INDEX[i]] == "native.pdb"):
                TOPTEN = 1
                NATIVEINDEX = i
                break
        if(NATIVEINDEX == None):
            NATIVEINDEX = len(ResultList) + 1
        f = open("./Results/" + TestDataset + "/Order.txt",'a')
        f.write(structurename + '\t' + str(NATIVEINDEX) + '\t' + str(pea) + '\n')
        f.close()

    if(options == 'S'):
        plt.figure()
        plt.scatter(DecoyRMSD,ResultList,color = "red")
        plt.xlabel("RMSD")
        plt.ylabel("Energy")
        plt.title(structurename)
        savefig("./Results/" + TestDataset + "/Pics/" + structurename + '.png')
        plt.close()
        EFile = open("./Results/" + TestDataset + "/Energy/Energy_" + structurename + ".txt","a")
        for i in range(len(ResultList)):
            EFile.write(str(ResultList[i]))
            EFile.write("\t")
            EFile.write(str(DecoyRMSD[i]))
            EFile.write("\n")
        EFile.close()
        pea = pearson(DecoyRMSD,ResultList)
        ENERGY = np.asarray(ResultList)
        INDEX = np.argsort(ENERGY)
        NATIVEINDEX = None
        for i in range(10):
            if(DecoyList[INDEX[i]] == "native.pdb"):
                TOPTEN = 1
                NATIVEINDEX = i
                break
        if(NATIVEINDEX == None):
            NATIVEINDEX = len(ResultList) + 1
        f = open("./Results/" + TestDataset + "/Order.txt",'a')
        f.write(structurename + '\t' + str(NATIVEINDEX) + '\t' + str(pea) + '\n')
        f.close()
