
import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.pyplot import savefig



def asis(x,y,matrix):
    return matrix[x][y]

def contourmap(amino1,amino2,layer):

    f = open("./predict/AAdistri.npy")

    cdDict={"ALA":{},"VAL":{},"LEU":{},"ILE":{},"PHE":{},\
            "TRP":{},"MET":{},"PRO":{},"GLY":{},"SER":{},\
            "THR":{},"CYS":{},"TYR":{},"ASN":{},"GLN":{},\
            "HIS":{},"LYS":{},"ARG":{},"ASP":{},"GLU":{},}

    layermatrix = [0,6,7,8,9,10]
    AAlist = cdDict.keys()
    AAlist.sort()
    for a1 in AAlist:
        for a2 in AAlist:
            cdDict[a1][a2] = np.load(f)

    matrix = cdDict[amino1][amino2][layer,:,:]
    #print matrix
    #print np.shape(matrix)
    x = y = range(0,20)
    #theta = range(0,180,9)
    #phi = range(0,360,18)
    #plt.figure()
    fig1 = plt.contourf(x,y,asis(x,y,matrix))
    #plt.grid(True)
    plt.colorbar(fig1)
    plt.xticks(x,fontsize=12)
    plt.yticks(y,fontsize=12)
    plt.title(amino1 + '-' + amino2 + "(" + str(layermatrix[layer]) + '-' + str(layermatrix[layer+1]) + r"$\AA$"+')', fontsize = 20)
    plt.xlabel(r'$\theta$',fontsize=14)
    plt.ylabel(r"$\phi$",fontsize=14)
    #plt.show()
    savefig("./predict/static/pics/" + amino1 + '-' + amino2 + '-' + str(layer) + ".png")
    #plt.close()

    


if __name__ == "__main__":
    args = sys.argv[1:]
    amino1 = args[0]
    amino2 = args[1]
    layer = int(args[2])
    contourmap(amino1,amino2,layer)
    

