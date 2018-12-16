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
from progress.bar import Bar
#import pylab

DO = 8.0 #distance cut-off, there is not any interaction below 8A
DELTA = 1.5 #parameter of the logistic probability function


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

def readChainPDB(filename):
    '''
    Gets the chains of the pdb
    '''
    try :
        f=open(filename,"r")
    except OSError :
        sys.exit("The file does not exist in the directory, please provide an existing file\n")
    else :
        list_chains = []
        line = f.readline()
        while not re.search("^COMPND", line) :
            line = f.readline()
        #The line contains COMPND information
        while re.search("^COMPND", line) :
            if re.search("CHAIN:", line) :
                list_chains.append(line.split()[-1][0])
            line = f.readline()
        f.close()
        return(list_chains)


def readPDB(filename, chain):
    f=open(filename,"r")
    atomes=[]
    re_end_chain=re.compile("^TER")
    #initialization of the first line
    line = f.readline()
    #The function searches lines with atoms
    while not (re.search("^ATOM",line) or re.search("^HETATM", line)) :
        line = f.readline()
    #The line is an atom, now the function searches for the right chain
    while str(line[21]) != chain :
        line = f.readline()
    #The line contains the right chain
    while not re_end_chain.search(line) : #If this is true, the chain ends
        if re.search("^ATOM",line) or re.search("^HETATM", line):
            if re.search("CA",line) and int(line[22:26])>0 :
                atomes.append(Atome(line))
        line = f.readline()
        
    f.close()
    return atomes

def dssp(filename, chain):
    '''
    Creates a list of secondary structures assignment 
    It reads a file out of DSSP and gets the 
    secondary structures
    '''
    list_ss = [] #list of secondary structures
    with open(filename, 'r') as dssp_file :
        line = dssp_file.readline()
        while not re.search('  #  RESIDUE', line) :
            #Reads until it reads a line with secondary structure
            line = dssp_file.readline()
        for line in dssp_file :
            if line[13]!= '!' and line[11] == chain: 
                #Adds the secondary structure
                list_ss.append((line[16]))
    return(list_ss)

def distance(atom1, atom2):
    '''
    Calculates and returns the distance between two Atome instances based on coordinates, using sqrt function from math module
    '''
    d=sqrt((atom1.xpos-atom2.xpos)**2+(atom1.ypos-atom2.ypos)**2+(atom1.zpos-atom2.zpos)**2)
    return d

def contacts_matrix(filename, DO, DELTA, chain) :
    '''
    Creates the contact matrix of the residus within the protein
    Returns an array of the size of the number of residues with probabilities
    of contacts
    '''
    list_atoms = readPDB(filename, chain) #the list of atoms from the pdb file
    contacts = np.zeros((len(list_atoms), len(list_atoms))) #initializes a matrix
    #of zeros the size of number of residues
    for i in range(len(list_atoms)) :
        for j in range(i,len(list_atoms)) : #the matrix is symetric
            #i is the line, j the col
            d = distance(list_atoms[i], list_atoms[j])
            contacts[i,j] = 1/((1+exp((d - DO) / DELTA)))
    return(contacts)

def single_PI(contacts, a, b, list_ss) :
    '''
    Calculates the PI cutting the contacts matrix between a and b
    If a or b cuts inside a secondary structure, it returns 0
    '''
    ss_a = list_ss[a]
    ss_b = list_ss[b]
    #If those are coils, it puts NA instead to facilitate the rest of the function
    if ss_a == " ":
        ss_a = "NA"
    if ss_b == " ":
        ss_b = "NA"
    if a == 0 and b != (len(list_ss)-1) and list_ss[b+1] == ss_b  :
        #The PU begins at the begining of the protein
        #The PU does not end at the end of the protein
        #Cutting at b cuts within a secondary structure
        return(0)
    elif a !=0 and b == (len(list_ss)-1) and list_ss[a-1] == ss_a :
        #The PU does not begin at the begining of the protein
        #The PU ends at the end of the protein
        #Cutting at a is within a secondary strucutre
        return(0)
    elif list_ss[b+1] == ss_b or list_ss[a-1] == ss_a :
        #The PU is inside the protein
        #Cutting at a or b cuts within a secondary structure
        return(0)
    else :
        A = 0;  B = 0;  C = 0
        for i in range(len(list_ss)) :
            for j in range(i, len(list_ss)) :
                if (i >= a and i <= b) and (j >= a and j<= b) : #the contact is in A
                    A += contacts[i,j]
                elif (i < a or i > b) and (j < a or j > b) : #the contact is in B
                    B += contacts[i,j]
                else : #the contact is between A and B
                    C += contacts[i,j]
        PI = (A * B - C**2)/((A + C) * (B + C))
        return(PI)

def calculate_PI(contacts, begin, min_size, max_size, list_ss) :
    '''
    Calculates PIs for a beginning of PU
    '''
    list_PI = []
    for end in range((begin+min_size), (begin+max_size)) :
        list_PI.append(single_PI(contacts, begin, end, list_ss))
    return(list_PI)





def main() :
    namefile = sys.argv[1]

    list_chains = readChainPDB(namefile)
    chain = input("Choose a chain, in your pdb file the chains are {} :"
        .format(','.join(list_chains)))
    while chain not in list_chains :
        print("It is not a chain from your pdb file")
        chain = input("Choose a chain, in your pdb file the chains are {} :"
        .format(','.join(list_chains)))
    name = (namefile.split("/")[-1]).split(".")[0]
    contacts = (contacts_matrix(namefile, DO, DELTA, chain))
    list_ss = dssp("DSSP/"+name+".out", chain)
    min_size = 9
    max_size = 39
    dico_PI = {}

    bar = Bar('Processing', max=(len(list_ss)-max_size))
    for begin in range(0,(len(list_ss)-max_size)) :
        dico_PI[begin+1] = calculate_PI(contacts, begin, min_size, max_size, list_ss)
        bar.next()
    bar.finish()
    print(dico_PI)
    #pylab.matshow((contacts))
    #pylab.colorbar()
    #pylab.show()


if __name__=='__main__':
    main()
