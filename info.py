#Notre but est de découper une protéine en PU et de trouver les meilleures zones de découpage.
#D'abord on va créer une matrice de contacts avec la fonction logistique : p(i,j)=1/((1+exp[d(i,j)−d0)/Δ])
#d(i,j) est la distance euclidienne entre deux carbones alpha
#d0 est la distance cut-off pour laquelle on considère qu'il n'y a plus de contact
#d0 est fixé à 8, Δ à 1.5

#Un PU c'est quand on a plus de contacts au sein de la PU que avec l'exterieur.

#Au début la protéine est entière : 1 grand PU compris entre i=1 et j=N
#Puis on découpe deux sous unités à 20 aa.
#On calcule le PI (PIi,j(m)=(AB−C^2)/((A+C)(B+C)))
#A est le nombre de contacts dans la sous unité A, B dans la sous unité B, et C entre A et B.
#m est la position entre i et j, celle qui va bouger. A la première itération est est de 20.

#Puis on passe à m=21, puis m=22, etc jusqu'à 60
#Puis on décale i de 1 et on recommence.

#Pour chaque i de départ, on aura un plot de PI en fonction de m. On définira un seuil.
#On pourra définir le meilleur PI.

#On va aussi utiliser un deuxième critère : de compacité ou ce séparation

#On lit un pdb, on récupère les coordonnées, on fait la matrice de contacts
#Avec les distances
#On decoupe en deux sous-unités, on teste et on calcule PI, compacité et séparation.
#On récupère les infos avec
#pos_i pos_j PI(m) compact(PUa, PUb) separation(PUa, PUb)


#On voudra le plus grand PI, le plus grand separation et le plus petit compact.

import re
import sys
from math import sqrt
from math import exp
import numpy as np

DO = 8.0 #distance cut-off
DELTA = 1.5


'''
class Atome : object with atoms informations
	the name of the atom
	the chain the atom comes from
	x,y and z coordinates of the atom
	residu type with 3 letters code name
	residu number in the protein
'''
class Atome :
	def __init__(self,line):
		self.atome_name=line[13:16].split(" ")[0]
		self.chain=str(line[21])
		self.xpos=float(line[30:38])
		self.ypos=float(line[38:46])
		self.zpos=float(line[46:54])
		self.residu_type=line[17:20]
		self.residu_num=int(line[22:26])
		self.atome_num=int(line[6:11])
	def __str__(self):
		return "{} {} {} {} {} {} {} {}".format(self.atome_name, self.chain, self.xpos, self.ypos , self.zpos, self.residu_type, self.residu_num, self.atome_num)

'''
readPDB(filename) fonction :
Reads the provided file in the data directory, uses only the first model
'''
def readPDB(filename):
	#files=os.listdir('../data/')
	try :
		f=open(filename,"r")
	except OSError :
		sys.exit("The file does not exist in the directory, please provide an existing file\n")
	else :
		atomes=[]
		re_atomes=re.compile("^ATOM")
		re_fin_modele=re.compile("^ENDMDL")
		for line in f:
			if re_atomes.search(line):
				if re.search("CA",line) :
					atomes.append(Atome(line))
			elif re_fin_modele.search(line):
				break
		f.close()
		return atomes

def distance(atom1, atom2):
    '''
    Calculates and returns the distance between two Atome instances based on coordinates, using sqrt function from math module
    '''
    d=sqrt((atom1.xpos-atom2.xpos)**2+(atom1.ypos-atom2.ypos)**2+(atom1.zpos-atom2.zpos)**2)
    return d

def contacts_matrix(filename, DO, DELTA) :
	'''
	Creates the contact matrix of the residus within the protein
	Returns an array of the size of the number of residues with probabilities
	of contacts
	'''
	list_atoms = readPDB(filename) #the list of atoms from the pdb file
	contacts = np.zeros((len(list_atoms), len(list_atoms))) #initializes a matrix
	#of zeros the size of number of residues
	for i in range(len(list_atoms)) :
		for j in range(len(list_atoms)) : #the matrix is symetric
			#i is the line, j the col
			d = distance(list_atoms[i], list_atoms[j])
			contacts[i,j] = 1/((1+exp((d - DO) / DELTA)))
	return(contacts)





def main() :
	namefile = sys.argv[1]
	print(contacts_matrix(namefile, DO, DELTA))




if __name__=='__main__':
    main()
