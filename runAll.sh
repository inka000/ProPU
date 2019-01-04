#!/bin/bash
display_help() {
    echo "Usage: $0 [option...]" >&2
    echo
    echo "   -h, --help                   Displays help"
    echo
    echo "   -i, --input                  Directory where pdb files are "
    echo
    echo "   --min                        Minimal size for a PU (10 by default)"
    echo
    echo "   --max                        Maximal size for a PU (40 by default)"
    echo
    echo "   --delta                      Parameter of the logistic probability function"
    echo "                                (1.5 by default)"
    echo
    echo "   --dist                       Distance cut-off for interactions"
    echo "                                (8.0 by default)"
    echo
    echo "   --all                        Whether it considers all PUs or only the best one"
    echo "                                [yes/no] (yes by default)"
    echo
    exit 1
}

run_program() {
    #####Repositories creation#####
    rm -Rf Query DSSP suivi
    mkdir Query DSSP suivi
    if ! [ -d "./resultPI" ]; then
        mkdir resultPI
    fi

    #####PDB queries will be relocalized in Query repository#####

    for query in "$queryFiles"* ; do
        cp $query ./Query
    done

    ######Secondary structure assignation with DSSP#####
    ######And criterions calculations ##################
    for f in Query/* ; do
        pdb=`echo $f | cut -d/ -f2`
        name=`echo ${f##*/} | cut -d. -f1`
        if ! [ -d "./resultPI/$name" ]; then
            mkdir resultPI/$name
        fi
        ./dssp-2.0.4-linux-amd64 -i ./$f -o DSSP/$name.out &>> suivi/tmp.log
        echo 'Query ='$name
        python3 calculationbis.py $f $dist $delta $minsize $maxsize $all
    done

    ######Analyses of PUs#######
    for f in resultPI/* ; do
        texte=`echo $f | cut -d/ -f2`
        name=`echo ${f##*/} | cut -d. -f1`
        echo 'Analyses of '$name
    done
}

while :
do
    case "$1" in
        -i | --input)
            if [ $# -ne 1 ] ; then 
                queryFiles="$2";
            else 
                break
            fi
            shift 2
            ;;
        --min)
            if [ $# -ne 1 ] ; then
                minsize="$2";
            else
                break
            fi
            shift 2
            ;;  
              
        --max)
            if [ $# -ne 1 ] ; then
                maxsize="$2";
            else
                break
            fi
            shift 2
            ;;
        --delta)
            if [ $# -ne 1 ] ; then
                delta="$2";
            else
                break
            fi
            shift 2
            ;;
        --dist)
            if [ $# -ne 1 ] ; then
                dist="$2";
            else
                break
            fi
            shift 2
            ;;
        --all)
            if [ $# -ne 1 ] ; then
                all="$2";
            else
                break
            fi
            shift 2
            ;; 
        -h | --help)
            display_help  # Call help function
            exit 0
            ;;
        --) # End of all options
            shift
            break
            ;;
        -*)
            echo "Error: Unknown option: $1" >&2
            display_help
            exit 1 
            ;;
        *)  # No more options
            break
            ;;
    esac
done

shift $((OPTIND - 1))


#Check parameters and set default values if parameters are empty
#or if parameters is not numeric
if [ -z ${minsize+x} ] || ! [ $(echo $minsize | grep -v [a-Z]) ]; then 
    minsize=10 
fi

if [ -z ${maxsize+x} ] || ! [ $(echo $maxsize | grep -v [a-Z]) ]; then 
    maxsize=40 
fi

if [ -z ${delta+x} ] || ! [ $(echo $delta | grep -v [a-Z]) ]; then 
    delta=1.5
fi

if [ -z ${dist+x} ] || ! [ $(echo $dist | grep -v [a-Z]) ]; then 
    dist=8.0 
fi

if [ -z ${all+x} ] ||  [ $(echo $dist | grep -v [a-Z]) ]; then 
    all="yes"
fi

echo "$minsize $maxsize $delta $dist"

#Check parameters values
if [ -z ${queryFiles+x} ]; then 
    echo "please provide a Directory"; 

elif [ $minsize -lt 0 ] || [ $maxsize -le 0 ] || [ $maxsize -lt $minsize ]; then
    echo "please provide correct min and max sizes"

elif [ $(bc <<< "$delta <= 0.0") -eq 1 ]; then
    echo "please provide a correct delta"

elif [ $(bc <<< "$dist <= 0.0") -eq 1  ]; then
    echo "please provide a correct distance"

else 
    run_program 
fi

