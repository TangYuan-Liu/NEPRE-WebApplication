import os
import numpy as np
import AminoAcid as AA
import math
import argparse

def Pearson(rmsd,energy):
    size = np.shape(rmsd)[0]
    x = np.empty(shape=[2,size])
    for i in range(size):
        x[0][i] = rmsd[i]
    for j in range(size):
        x[1][j] = energy[j]
    y = np.corrcoef(x)
    return y[0][1]


def load_EnergyMatrix(cutoff,base_path="./"):
    aaDict={"ALA":{},"VAL":{},"LEU":{},"ILE":{},"PHE":{},\
            "TRP":{},"MET":{},"PRO":{},"GLY":{},"SER":{},\
            "THR":{},"CYS":{},"TYR":{},"ASN":{},"GLN":{},\
            "HIS":{},"LYS":{},"ARG":{},"ASP":{},"GLU":{},}

    List = aaDict.keys()
    List.sort()
    
    f = open(base_path + str(cutoff) + ".npy")
    for amino1 in List:
        for amino2 in List:
            aaDict[amino1][amino2] = np.load(f)
    f.close()
    return aaDict

def extract_Data(line):
    """
    This part will extracted data from line according to the standard 
    PDB file format(Version 3.3.0, Nov.21, 2012)
    """
    res = []

    line = line.strip()
    #record_name
    res.append(line[0:4].strip(' '))

    #atom_serial
    res.append(line[6:11].strip(' '))

    #atom_name
    res.append(line[12:16].strip(' '))

    #alternate_indicator
    res.append(line[16])

    #residue_name
    res.append(line[17:20].strip(' '))

    #chain_id
    res.append(line[21].strip(' '))

    #residue_num
    res.append(line[22:26].strip(' '))

    #xcor
    res.append(line[30:38].strip(' '))

    #ycor
    res.append(line[38:46].strip(' '))

    #zcor
    res.append(line[46:54].strip(' '))
   
    return res
    

def calculate_Energy(f,matrix,cutoff):

    CurrentAANitrogen = None
    CurrentAACA = None
    Currentresidue_num = None
    EachAA = []    
    CurrentAA = None 

    for line in f.readlines():        
        if(line[0:4] != "ATOM"):
            continue
        element_list = extract_Data(line)
        record_name = element_list[0]
        atom_name = element_list[2]
        residue_name = element_list[4]
        alternate_indicator = element_list[3]
        residue_num = element_list[-4]
        xcor = float(element_list[-3])
        ycor = float(element_list[-2])
        zcor = float(element_list[-1])
        
        if(atom_name == "H"):
            continue
        if(residue_name not in matrix):
            continue

        if(CurrentAA == None):
            CurrentAA = AA.AminoAcid(residue_name)
            Currentresidue_num = residue_num
            if(atom_name == "N" or atom_name == "CA"):
                if(alternate_indicator == "B"):
                    continue
                if(atom_name == "N"):
                    CurrentAANitrogen = np.array([xcor,ycor,zcor])
                else:
                    CurrentAACA = np.array([xcor,ycor,zcor])
            if(residue_name == "GLY" or atom_name not in {"N","CA","C","O","O1","02"}):
                if(alternate_indicator != " "):
                    #If cases like "AASN or BASN" appears, we only add A 
                    if(alternate_indicator == "A" and line[15] == "1"):
                        CurrentAA.SumCenters(xcor,ycor,zcor)
                    else:
                        continue
                else:
                    CurrentAA.SumCenters(xcor,ycor,zcor)
        else:
            #If another amino acid begins
            if(residue_num != Currentresidue_num):
                state = CurrentAA.CalculateCenter()
                if(state == False):
                    CurrentAA = AA.AminoAcid(residue_name)
                    Currentresidue_num = residue_num
                    continue

                CurrentAA.InputCAN(CurrentAANitrogen,CurrentAACA)
                #Amino Acid check
                EachAA.append(CurrentAA)
                del CurrentAA
                CurrentAA = AA.AminoAcid(residue_name)

                Currentresidue_num = residue_num
                if(atom_name == "N" or atom_name == "CA"):
                    if(alternate_indicator == "B"):
                        continue
                    if(atom_name == "N"):
                        CurrentAANitrogen = np.array([xcor,ycor,zcor])
                    else:
                        CurrentAACA = np.array([xcor,ycor,zcor])
                if(residue_name == "GLY" or atom_name not in {"N","CA","C","O","O1","02"}):
                    if(alternate_indicator != " "):
                    #If cases like "AASN or BASN" appears, we only add A 
                        if(alternate_indicator == "A" and line[15] == "1"):
                            CurrentAA.SumCenters(xcor,ycor,zcor)
                        else:
                            continue
                    else:
                        CurrentAA.SumCenters(xcor,ycor,zcor)
            #If still the same amino acid
            else:
                if(atom_name == "N" or atom_name == "CA"):
                    if(alternate_indicator == "B"):
                        continue
                    if(atom_name == "N"):
                        CurrentAANitrogen = np.array([xcor,ycor,zcor])
                    else:
                        CurrentAACA = np.array([xcor,ycor,zcor])
                if(residue_name == "GLY" or atom_name not in {"N","CA","C","O","O1","02"}):
                    if(alternate_indicator != " "):
                    #If cases like "AASN or BASN" appears, we only add A 
                        if(alternate_indicator == "A" and line[15] == "1"):
                            CurrentAA.SumCenters(xcor,ycor,zcor)
                        else:
                            continue
                    else:
                        CurrentAA.SumCenters(xcor,ycor,zcor)
    
    state = CurrentAA.CalculateCenter()
    if(state != False):        
        CurrentAA.CalculateCenter()
        CurrentAA.InputCAN(CurrentAANitrogen,CurrentAACA)
        EachAA.append(CurrentAA)
    #Scan over. Each amino acid is stored as an object in EachAA. Next step is to calculate the energy, results will be saved in EnergyList. 
    
    #Store the energy
    E = 0 
    for m in range(len(EachAA)):
        #Establish axis first    
        EachAA[m].EstablishCoordinate()
        for n in range(len(EachAA)):
            if(m == n):
                continue
            else:
                dis = EachAA[m].DistanceBetweenAA(EachAA[n].center)
                
                if(dis < cutoff):#If the distance between two amino acid less than 10, we believe the two amino acid have interaction
                    rho,theta,phi = EachAA[m].ChangeCoordinate(EachAA[n].center)
                    theta = min(int(math.floor(theta*20/np.pi)),19)
                    phi = min(int(math.floor(phi*10/np.pi) + 10),19)
                    E += matrix[EachAA[m].name][EachAA[n].name][theta][phi] / rho 
                 
    return E


if __name__ == "__main__":
     
    parser = argparse.ArgumentParser(description="Nepre-F Scoring Function Created by CSRC")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-s","--single",help="calculate single PDB",action="store_true")
    group.add_argument("-m","--multi",help="calculate a series of PDB",action="store_true")
    parser.add_argument("-o","--output",help="save the results as a text file in running folder",action="store_true")
    parser.add_argument("path",help="PDB file path of folder path")
    parser.add_argument("cutoff",type=int,help="cutoff parameter for Nepre-F")
    args = parser.parse_args()
    
    if(args.single == True):
        c = args.cutoff
        matrix = load_EnergyMatrix(c)
        p = args.path
        f = open(p)
        E = calculate_Energy(f,matrix,c)
        print "Nepre Potential Energy"
        print "Using Cutoff:",c
        print p,E
        if(args.output):
            save_file = open("./latest_results.txt","wb")
            save_file.write("Nepre Potential Energy" + '\n')
            save_file.write("Using Cutoff:" + str(c) + '\n')
            save_file.write(p + '\t' + str(E))
            save_file.close()
    if(args.multi == True):
        c = args.cutoff
        matrix = load_EnergyMatrix(c)
        folder_path = args.path
        file_list = []
        for pdb_file in os.listdir(folder_path):
            file_list.append(pdb_file)
        E = []
        if(folder_path[-1] != '/'):
            folder_path += '/'
        for pdb_file in file_list:
            pdb_path = folder_path + pdb_file
            f = open(pdb_path)
            E.append(calculate_Energy(f,matrix,c))
        if(args.output):
            save_file = open("./latest_results.txt","wb")
            save_file.write("Nepre Potential Energy" + '\n')
            save_file.write("Using Cutoff:" + str(c) + '\n')
            for i in range(len(E)):
                save_file.write(file_list[i] + '\t' + str(E[i]))
                save_file.write('\n')
            save_file.close()
    
        print "Nepre Potential Energy"
        print "Using Cutoff:",c
        for i in range(len(E)):
            print file_list[i],'\t',E[i]
