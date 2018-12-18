import matplotlib.pyplot as plt
import sys


def read_info(namefile) :
    '''
    Reads the file containing information about PUs
    Returns three lists :
    - list of PI
    - list of sigma
    - list of k
    '''
    list_PI = []
    list_sigma = []
    list_k = []
    list_size = []
    list_begin = []
    with open(namefile, "r") as file :
        line = file.readline()
        for line in file :
            list_begin.append(float(line.split("\t")[0]))
            list_size.append(float(line.split("\t")[1]))
            list_PI.append(float(line.split("\t")[2]))
            list_sigma.append(float(line.split("\t")[3]))
            list_k.append(float(line.split("\t")[4]))
    return(list_begin, list_size, list_PI, list_sigma, list_k)







def main() :
    namefile = sys.argv[1]
    read = read_info(namefile)
    list_begin = read[0]
    list_size = read[1]
    list_PI = read[2]
    list_sigma = read[3]
    list_k = read[4]

    list_temp_size = []
    list_temp_PI = []
    for i in range(len(list_begin)) :
        if list_begin[i] == "277" :
            print("ok")
            list_temp_PI.append(float(list_PI[i]))
            list_temp_size.append(float(list_size[i]))


    #plt.plot(list_temp_size, list_temp_PI, 'ro')
    #plt.axis([min(list_temp_size), max(list_temp_size), min(list_temp_PI), max(list_temp_PI)])
    #plt.show()

    plt.plot(list_size, list_k, 'ro')
    plt.axis([min(list_size), max(list_size), min(list_k), max(list_k)])
    plt.show()


if __name__=='__main__':
    main()