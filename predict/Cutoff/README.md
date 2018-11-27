# Nepre-Potential(Cutoff)
Scoring Function based on Neighbourhood Preference Statistics  

Usage
----------
The runing folder should contain:
* AminoAcid.py(Class for establish amino acid)
* cutoff.npy(Energy matrix)

We provide 7 cutoff options with cutoff between 4 angstrom and 10 angstrom.

For single protein potential energy calculate, choose a cutoff(6 angstrom e.g) and go to linux shell and type:
<pre><code>
Nepre@liulab:~$ python Nepre.py 6 ./example.pdb
</code></pre>

For multi-object job, you could speed up by using Nepre as a module:

<pre></code>
import Nepre

#choose a cutoff
cutoff = 6

#select a protein
path = "./example.pdb"
f = open(path)

#load energy matrix
Matrix = Nepre.load_EnergyMatrix(cutoff)

#calculate Nepre potential energy
E = Nepre.calculate_Energy(f,Matrix,cutoff)
</code></pre>

Extensions
----------
Nepre module also provide some useful function:
* Calculate the pearson coefficient correlation.
* Extract data from standard PDB file.
<pre><code>
"""
Pearson Coefficient
"""
import Nepre
x = [1,2,3,4]
y = [1,2,3,4]
p = Nepre.Pearson(x,y)

"""
Extract Data
"""
import Nepre
f = open("./example.pdb")
res = []
for line in f.readlines():
    res.append(Nepre.extract_Data(line))
</code></pre>
