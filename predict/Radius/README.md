# Nepre-Potential(Radius)
Scoring Function based on Neighbourhood Preference Statistics  

Usage
----------
The runing folder should contain:
* AminoAcid.py(Class for establish amino acid)
* radius.npy(Energy matrix)

A statistic of radius of different kind of amino acid have been done by us. We use gaussian distribute function to 
fit them and use the gaussian mean data as the default vlaue. See details in mean_radius.txt.

For single protein potential energy calculate, go to linux shell and type:
<pre><code>
Nepre@liulab:~$ python Nepre.py ./example.pdb
</code></pre>

For multi-object job, you could speed up by using Nepre as a module:

<pre></code>
import Nepre

#select a protein
path = "./example.pdb"
f = open(path)

#load energy matrix
Matrix = Nepre.load_EnergyMatrix()

#calculate Nepre potential energy
E = Nepre.calculate_Energy(f,Matrix)
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
