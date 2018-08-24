#/home/jjuarez/soft/iEBDD/accessory/miniconda3/bin//python

###     CleanPDB v 1.0    ###
#  University of Edinburgh  #
#        July 2018          #
#############################

__author__ = "Jordi Juarez-Jimenez"


errorcode = 999

# Imports
import sys, os


# Subroutines #



def clean_pdb(pdbfile):
    '''

    :param pdbfile:(file opened with reading permissions)
    :return: (int)

    Rewrites pdbfile ignoring deleting all entries different than "ATOM" "TER" and "END"

    '''

    try:
        with open('_tempfile.pdb', 'w') as f:

            while True:
                line = pdbfile.readline()
                if line == '':
                    break
                else:
                    tmpline = line.replace("1A", "1 ")
                    line = tmpline.replace("2A", "2 ")
                    mor = line.split()
                    if mor[0] == "ATOM" or mor[0] == "TER" or mor[0] == "END":
                        if "1B" not in line and "2B" not in line:
                            f.write(line)


        return 0

    except IOError:
        message = "cleanPDB error: It was not possible to complete rewriting the pdb file"
        print(message)
        return errorcode



def split_chains(pdbfile, outrootname):
    '''

    :param pdbfile: (file opened with reading permissions)
    :return: (int)

    Rewrites pdbfile splitting it in one file for each chain

    '''

    allchains = {}
    currentchain = ""
    try :
        while True:
            line = pdbfile.readline()
            if line == '':
                break
            else:
                mor = line.split()
                if len(mor) > 4 and mor[0] != "TER" and mor[4] not in allchains:
                    allchains[mor[4]] = []
                    currentchain = mor[4]
                allchains[currentchain].append(line)

        for chain in allchains:
            with open(outrootname + "_%s.pdb" % chain , 'w') as f:
                for line in allchains[chain]:
                    f.write(line)

        return 0

    except IOError:
        raise
        return errorcode





if __name__ == "__main__":

    if len(sys.argv) != 2:
        print ("Usage: cleanPDB.py pdbfile")
        sys.exit()


    filename = sys.argv[1]
    message = "Procesing %s..." % filename
    print(message)

    outrootname = filename.split('.')[0]

    with open(filename) as inpf:
        status = clean_pdb(inpf)
        if status != 0:
            message = "Unfortunately cleanPDB failed to process the following pdb file: %s" % filename
            print(message)
            sys.exit(status)
        else:
            with open('_tempfile.pdb') as tmpfile:
                status = split_chains(tmpfile, outrootname)

            command="rm _tempfile.pdb %s" %filename
            os.system(command)
            message="Done!"
            print(message)
            sys.exit(status)
