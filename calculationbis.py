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
from progress.bar import FillingSquaresBar
import matplotlib.pyplot as plt
import random
import statistics
from scipy.stats import norm




class Atome :
    '''
    class Atome : This class groups atoms informations

    Attributes :
        atome_name (str) : the name of the atom
        chain (str) : the chain the atom comes from
        xpos, ypox and zpos (float) : x,y and z coordinates of the atom
        residu_type (str) : residu type with 3 letters code name
        residu_num (int) : residu position in the protein
        atome_num (int) : atom position in the protein
    '''
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
        '''
        Returns formated expression of the atom informations
        '''
        print(" name : {} \n chain : {} \n xpos :{} ypos : {} zpos : {} \nresidu_type : {} \nresidu_num : {} \natome_num : {}".format(self.atome_name, self.chain, self.xpos, self.ypos , self.zpos, self.residu_type, self.residu_num, self.atome_num))
    

    def distance(self, atom2) :
        '''
        Calculates and returns the distance between two Atome instances based on coordinates, using sqrt function from math module
        '''
        d=sqrt((self.xpos-atom2.xpos)**2+(self.ypos-atom2.ypos)**2+(self.zpos-atom2.zpos)**2)
        return(d)


class PU :
    '''
    class PU :  This class groups Protein Unit informations

    Attributes :
        begin (int) : position of the begining of the PU
        size (int) : size in aa of the PU
        PI (float) : partition index of the PU
        Sigma (float) : separation criterion
        k (float) : compactness criterion
    '''
    def __init__(self, begin, end, size) :
        self.begin = begin
        self.end = end
        self.size = size
        self.PI = None
        self.sigma = None
        self.k = None
        self.signif = None

    def __str__(self):
        '''
        Returns formated expression of the PU informations
        '''
        print("begin : {} \nend : {} \nsize : {} \nPI : {} \nsigma : {} \nk : {} \nsignif : {}\n".format(self.begin, self.end, self.size, self.PI , self.sigma, self.k, self.signif))

    def single_Sigma(self, contacts, a, b, flag) :
        '''
        Calculates the separation criterion for the PU between a and b.
        A low sigma means that the PU does less contacts with the rest of the protein
        '''
        alpha = 0.43 #from the article
        Pinter = 0 #interactions of A vs the rest of the protein
        Ptot = 0 #total of interactions
        if flag == 0 : #The PI was not calculated
            return (0) 
        for i in range(contacts.shape[0]) :
                for j in range(i, contacts.shape[1]) :
                    if (i >= a and i <= b) and (j >= a and j<= b) : #the contact is in A
                        Ptot += contacts[i,j]
                    elif (i < a or i > b) and (j < a or j > b) : #the contact is in B
                        Ptot += contacts[i,j]
                    else : #the contact is between A and B
                        Ptot += contacts[i,j]
                        Pinter += contacts[i,j]
        haut = Pinter /( (self.size ** alpha) * (contacts.shape[0] - self.size)**alpha )
        bas = Ptot / contacts.shape[0]
        sigma = haut / bas
        self.sigma = sigma


    def single_criterion(self, contacts, a, b, list_ss) :
        '''
        Calculates the PI cutting the contacts matrix between a and b
        If a or b cuts inside a secondary structure, it returns 0
        Returns also the sigma (separation criterion) and the k 
        (compactness criterion)
        '''
        flag = 0
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
            return(0,0,0)
        elif a !=0 and b == (len(list_ss)-1) and list_ss[a-1] == ss_a :
            #The PU does not begin at the begining of the protein
            #The PU ends at the end of the protein
            #Cutting at a is within a secondary strucutre
            return(0,0,0)
        elif list_ss[b+1] == ss_b or list_ss[a-1] == ss_a :
            #The PU is inside the protein
            #Cutting at a or b cuts within a secondary structure
            return(0,0,0)
        else :
            flag = 1 #a PI is calculated
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
            self.single_Sigma(contacts, a, b, flag)
            k = A / (self.size)
            self.PI = PI
            self.k = k

    def add_signif(self, criterion) :
        '''
        Adds the criterion(s) that is/are significant
        - No letter means it is not a PU
        - P means that the PI is significantly high
        - S mean that the sigma is significantly low
        - K mean that k is significantly high
        '''
        self.signif = criterion




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
                l = line.split(":")[-1].split()
                for c in l :
                    list_chains.append(c[0])
            line = f.readline()
        f.close()
        return(list_chains)


def readPDB(filename, chain):
    '''
    Reads a pdb and create a list of alpha carbon atoms (Atome instances)
    '''
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
                ss = line[16]
                if ss == "T" or ss == "S" : 
                    #coils
                    ss = " "
                elif ss == "H" or ss == "G" or ss == "I" :
                    #helix
                    ss = "H"
                elif ss == "E" or ss == "B":
                    #b sheets
                    ss = "E"
                list_ss.append(ss)
    return(list_ss)


def contacts_matrix(filename, DO, DELTA, chain, list_atoms) :
    '''
    Creates the contact matrix of the residus within the protein
    Returns an array of the size of the number of residues with probabilities
    of contacts
    '''
    contacts = np.zeros((len(list_atoms), len(list_atoms))) #initializes a matrix
    #of zeros the size of number of residues
    for i in range(len(list_atoms)) :
        for j in range(i,len(list_atoms)) : #the matrix is symetric
            #i is the line, j the col
            d = list_atoms[i].distance(list_atoms[j])
            contacts[i,j] = 1/((1+exp((d - DO) / DELTA)))
    return(contacts)




def calculate_criterions(contacts, begin, min_size, max_size, list_ss) :
    '''
    Calculates PIs for a beginning of PU
    '''
    list_PU = []
    for end in range((begin+min_size), (begin+max_size)) :
        size = end - begin + 1
        pu = PU(begin+1, end+1, size)
        pu.single_criterion(contacts, begin, end, list_ss)
        if pu.PI != None :
            list_PU.append(pu)

    return(list_PU)


def calculate_zscore(PUs, stype) :
    '''
    Calculates the zscores for a list of given values of specified type
    It allows to normalize values to fit a normal distribution
    '''
    if stype == "PI" :
        list_score = [x.PI for x in PUs]
    elif stype == "sigma" :
        list_score = [x.sigma for x in PUs]
    else :
        list_score = [x.k for x in PUs]

    mscore = statistics.mean(list_score) #mean
    sdscore = statistics.stdev(list_score)  #standard deviation

    list_zscore = [(x-mscore)/sdscore for x in list_score]
    return(list_zscore)

def calculate_pvalue(list_score) :
    '''
    Calculates the p-values for a given list of zscores
    '''
    list_pvalues = []
    for score in list_score :
        if score >= 0 :
            list_pvalues.append(2*norm.cdf(-score)) #with normal distribution
        else :
            list_pvalues.append(2*norm.cdf(score))
    return(list_pvalues)



def find_PU(PUs, option) :
    '''
    #If option is not :
        Finds PUs based on PI and sigma :
        - the p-value of PI should be lower than 0.05 and the z-score should be higher than 0
        - the p-value of sigma should be lower than 0.05 and the z-score should be lower than 0
        - k is informative. The highest is the k, the better
    #If option is yes :
        Adds significativity based on PI, sigma and k values but does not sort PUs
    '''

    #Calculates z-scores
    list_PI = calculate_zscore(PUs, "PI")
    list_sigma = calculate_zscore(PUs, "sigma")
    list_k = calculate_zscore(PUs, "k")

    #Calculates p-values
    list_PI_pvalues = calculate_pvalue(list_PI)
    list_sigma_pvalues = calculate_pvalue(list_sigma)
    list_k_pvalues = calculate_pvalue(list_k)

    seuil = 0.05

    cpt_PU = 0
    best_PU = []
    for pval in list_PI_pvalues :
        #if pval < seuil and list_PI[cpt_PU] > 0 : 
        if pval < seuil and list_sigma_pvalues[cpt_PU] < seuil and list_k_pvalues[cpt_PU] < seuil :
            if list_PI[cpt_PU] > 0 and list_sigma[cpt_PU] < 0 and list_k[cpt_PU] > 0 :
                #PI, sigma and k are significant and none of PI, sigma and k are too low or high
                PUs[cpt_PU].add_signif("PSK")
        elif (pval < seuil and list_sigma_pvalues[cpt_PU] < seuil) :
            if list_PI[cpt_PU] > 0 and list_sigma[cpt_PU] < 0 :
                PUs[cpt_PU].add_signif("PS")
        elif (pval < seuil and list_k_pvalues[cpt_PU] < seuil) :
            if list_PI[cpt_PU] > 0 and list_k[cpt_PU] > 0 :
                PUs[cpt_PU].add_signif("PK")
        elif pval < seuil and list_PI[cpt_PU] > 0 :
            PUs[cpt_PU].add_signif("P")
        if option == "yes" :
            #it takes all PUs
            best_PU.append(PUs[cpt_PU])
        else :
            #It takes only the significant ones
            if PUs[cpt_PU].signif != None :
                best_PU.append(PUs[cpt_PU])
        cpt_PU += 1
    return(best_PU)

def single_best_PU(list_to_analyse) :
    '''
    Analyses a list of PU and returns the best one
    The best PU is the one with :
    - all significant and best criterions
    or
    - two significant and best criterions
    or
    - best PI
    '''
    list_PSK = []
    list_bi = []
    list_P = []
    list_none = []
    for index, pu in enumerate(list_to_analyse) :
        if pu.signif == "PSK" :
            list_PSK.append([pu,index])
        elif pu.signif == "PK" or pu.signif == "PS" :
            list_bi.append([pu,index])
        elif pu.signif == "P" :
            list_P.append([pu,index])
        else :
            list_none.append([pu, index])
    
    if len(list_PSK) == 1 :
        #There is just one PU with all criterions
        return (list_PSK[0][0], list_PSK[0][1])

    elif len(list_PSK) > 0:
        #There are several PUs with all criterions
        #The best PI is kept
        max_PI = list_PSK[0]
        for element in list_PSK :
            if element[0].PI > max_PI[0].PI :
                max_PI = element
        return(max_PI[0], max_PI[1])

    elif len(list_bi) == 1 :
        #there is only one PU with two criterions
        return(list_bi[0][0], list_bi[0][1])
    elif len(list_bi) > 0:
        #there are several PUs with two criterions
        #the best PI is kept
        max_PI = list_bi[0]
        for element in list_bi :
            if element[0].PI > max_PI[0].PI :
                max_PI = element
        return(max_PI[0], max_PI[1])

    elif len(list_P) == 1 :
        #There is only one PU with significant PI
        return(list_P[0][0], list_P[0][1])
    elif len(list_P) > 0 :
        #There are several PUs with significant PI
        #The best PI is kept
        max_PI = list_P[0]
        for element in list_P :
            if element[0].PI > max_PI[0].PI :
                max_PI = element
        return(max_PI[0], max_PI[1])

    elif len(list_none) == 1 :
        #there is not any significant criterion
        return(list_none[0][0], list_none[0][1])
    else :
        max_PI = list_none[0]
        for element in list_none :
            if element[0].PI > max_PI[0].PI :
                max_PI = element
        return(max_PI[0], max_PI[1])


def best_PU(PUs) :
    '''
    Finds best PUs based on criterions
    The functions calls single_best_PU to find the best PU
    of a list of PUs to analyse.
    To create a list of PUs to analyse :
    - it gets the first best PU
    - it gathers remaining PU that do not overlap the best PU
    - it gets the second best PU
    - it continues the same way until the list to analyse is empty
    '''
    #First, it gets the best PU of the protein
    best, index = single_best_PU(PUs)
    list_best_PU = [best]
    #Then, it creates a list of PUs to analyse 
    #that contains PU outside the first one
    list_to_analyse = [pu for pu in PUs if (pu.end < best.begin or pu.begin > best.end)]
    
    while len(list_to_analyse) > 0 :
        #there are other PU outside
        best, index = single_best_PU(list_to_analyse)
        list_best_PU.append(best)
        #The new list to analyse will contain PU outside the last best one
        list_to_analyse = [pu for pu in list_to_analyse if pu.end < best.begin or pu.begin > best.end]  

    return(list_best_PU)




def create_file(list_PU, name, namefile, min_size, option, chain) :
    '''
    Creates the file containing the calculated values :
    - PI
    - sigma
    - k
    '''
    file = open("./resultPI/"+name+"/"+namefile, option)
    file.write("Results for the chain {} of the protein {}\n".format(chain, name))
    file.write("begin\tend\tsize\tPI\tsigma\tk\tsignificant\n")
    for pu in list_PU :
        file.write("{:<4}\t{:<4}\t{:<4}\t{:<5}\t{:<5}\t{:<5}\t{:<3}\n".format(pu.begin, pu.end, \
            pu.size, round(pu.PI,3), round(pu.sigma,3), round(pu.k,3), str(pu.signif)))
    file.close()


def main() :
    namefile = sys.argv[1]
    name = (namefile.split("/")[-1]).split(".")[0]

    #parameters
    DO = float(sys.argv[2]) #distance cut-off, there is not any interaction below 8A
    DELTA = float(sys.argv[3]) #parameter of the logistic probability function
    MIN_SIZE = int(sys.argv[4]) #minimal size of a PU
    MAX_SIZE = int(sys.argv[5]) #maximal size of a PU
    option = str(sys.argv[6])#whether it considers all PUs or only significant ones

    #First, the programm reads the available chains from the pdb and ask the user
    #which chain he/she wants to analyse
    list_chains = readChainPDB(namefile)
    chain = input("Choose a chain, in your pdb file the chains are {} :"
        .format(','.join(list_chains)))
    while chain not in list_chains :
        print("It is not a chain from your pdb file")
        chain = input("Choose a chain, in your pdb file the chains are {} :"
        .format(','.join(list_chains)))


    list_atoms = readPDB(namefile, chain) #the list of atoms from the pdb file
    #Then, the programm creates the contacts matrix 
    contacts = contacts_matrix(namefile, DO, DELTA, chain, list_atoms)
    #It assigns secondary structures
    list_ss = dssp("DSSP/"+name+".out", chain)

    #And calculates PI, sigma and k criterions
    list_PU = []
   
    bar = FillingSquaresBar('Processing', max=(len(list_ss)-MAX_SIZE))
    for begin in range(0,(len(list_ss)-MAX_SIZE)) :
        list_PU = list_PU + calculate_criterions(contacts, begin, MIN_SIZE, MAX_SIZE, list_ss)
        bar.next()
    bar.finish()

    #Finds the best PU based on criterions values
    found_PU = find_PU(list_PU, option)

    create_file(found_PU, name, chain+"_"+name+"2.txt", MIN_SIZE, "w", chain)

    #best_PU(found_PU, MAX_SIZE)
    list_best_PU = best_PU(found_PU)

    create_file(list_best_PU, name, chain+"_"+name+".txt", MIN_SIZE, "w", chain)
    #for pu in list_PU :
    #    print(pu.__str__())
    #create_file(dico_PI, name+".txt", MIN_SIZE, "w")

    for pu in list_best_PU :
        pu.__str__()
        #plt.imshow(contacts)
        #plt.colorbar()
        #plt.axvline(pu.begin)
        #plt.axvline(pu.end)
        #plt.axhline(pu.begin)
        #plt.axhline(pu.end)
        #plt.text(pu.begin+2, pu.end-2, "A", fontsize = 20, color = "red")
    
        #plt.savefig("resultPI/"+name+"/"+chain+"_"+name+"_"+str(pu.begin)+".png")
        #plt.show()

    


if __name__=='__main__':
    main()
