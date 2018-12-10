from info import *

#Test of the pdb reading
def func1(filename) :
    return(readPDB(filename))

def test1() :
    atomes = func1("1BGXT.pdb")
    for at in atomes :
        assert at.atome_name == "CA"

#Test of the distance calculation
def func2(filename) :
    atomes = readPDB(filename)
    list_dist = []
    for i in range(len(atomes)) :
        for j in range(i, len(atomes)) :
            list_dist.append(distance(atomes[i], atomes[j]))
    return(list_dist)

def test2() :
    list_dist = func2("1BGXT.pdb")
    for dist in list_dist :
        assert dist >= 0

#Test of the contacts matrix size
def func3(filename) :
    return(contacts_matrix(filename, DO, DELTA).shape)

def test3() :
    assert func3("1BGXT.pdb")[0] == 828
    assert func3("1BGXT.pdb")[1] == 828

#Test of the diagonal of the contacts matrix
def func4(filename) :
    return(contacts_matrix(filename, DO, DELTA))

def test4() :
    mat = func4("1BGXT.pdb")
    for i in range(mat.shape[0]) :
        assert mat[i,i] == 1/((1+exp((0.0 - DO) / DELTA)))
