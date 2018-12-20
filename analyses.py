import matplotlib.pyplot as plt
import sys



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
    def __init__(self, line) :
        self.begin = float(line.split("\t")[0])
        self.size = float(line.split("\t")[1])
        self.PI = float(line.split("\t")[2])
        self.sigma = float(line.split("\t")[3])
        self.k = float(line.split("\t")[4])

    def __str__(self):
        '''
        Returns formated expression of the PU informations
        '''
        print("begin : {} \nsize : {} \nPI : {} \nsigma : {} \nk : {} ".format(self.begin, self.size, self.PI , self.sigma, self.k))
    

def read_info(namefile) :
    '''
    Reads the file containing information about PUs
    Creates a list of PUs (PU instances)
    '''
    PUs = []
    with open(namefile, "r") as file :
        line = file.readline()
        for line in file :
            PUs.append(PU(line))
    return(PUs)


def find_PU(PUs) :
    '''
    Finds the best PUs based on PI, sigma and k values :
    - PI should be in the interval [max(PI)-(0.1 * max(PI)), max(PI)]
    - sigma should be in the interval [min(sigma) , min(sigma) + (0.1 * min(sigma))]
    - k should be in the interval [max(k) - (0.1 * max(k)) , max(k)]
    '''
    #Gets the informations in lists
    #list_size = [x.size for x in PUs]
    #list_begin = [x.begin for x in PUs]

    #Finds best PI, sigma and k values
    best_PI = max([x.PI for x in PUs])
    print(best_PI)
    best_sigma = min([x.sigma for x in PUs])
    print(best_sigma)
    best_k = max([x.k for x in PUs])
    print(best_k)

    print([x for x in PUs if x.PI == best_PI][0].__str__())
    print("\n")
    print([x for x in PUs if x.sigma == best_sigma][0].__str__())
    print("\n")
    print([x for x in PUs if x.k == best_k][0].__str__())

    #Finds best PUs
    best_PU = []
    for pu in PUs :
        if (pu.k >= best_k - (best_k * 0.1)) :
            best_PU.append(pu)
            print("ok")

    return(best_PU)




def main() :
    namefile = sys.argv[1]
    PUs = read_info(namefile)
    


    #plt.plot(list_temp_size, list_temp_PI, 'ro')
    #plt.axis([min(list_temp_size), max(list_temp_size), min(list_temp_PI), max(list_temp_PI)])
    #plt.show()

    best_PU = find_PU(PUs)
    for pu in best_PU :
        print(pu.__str__())

    #plt.plot(list_begin, list_PI)
    #plt.axis([min(list_begin)-1, max(list_begin)+1, 0, 1])
    #plt.show()

    #plt.hist(list_PI)
    #plt.show()

if __name__=='__main__':
    main()