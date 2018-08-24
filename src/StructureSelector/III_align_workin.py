#!/usr/bin/python

import os, sys
import pymol, pickle

def calc_rmsd(mobpdb, refpdb, align_str):
    '''

    :param mobpdb: (str) A pdb file name (mobile structure)
    :param refpdb: (str) A pdb file name (reference structure)
    :param align_str: (str) A pymol like selection mask.
    :return: rmsd: (float) Not fitted rmsd value

    Uses pymol functions to calculate the rmsd among atoms in align_str between the mobpdb and refpdb

    '''

    pymol.pymol_argv = ['pymol', '-qc'] + sys.argv[1:]
    pymol.finish_launching()
    cmd = pymol.cmd


    cmd.load(mobpdb)
    cmd.load(refpdb)

    mobile_obj = mobpdb.rstrip('.pdb')
    ref_obj = refpdb.rstrip('.pdb')


    ca_rmsd = cmd.align("polymer and name CA and (%s)" % mobile_obj, "polymer and name CA and (%s)" % ref_obj,
                        object="test4", reset=0)

    rmsd = cmd.align("resi %s and (%s)" % (align_str, mobile_obj),
                     "resi %s and (%s)" % (align_str, ref_obj),  transform =0, object="test4", reset=1)

    cmd.delete(mobile_obj)
    cmd.delete(ref_obj)

    #print(ca_rmsd)
    #print(rmsd)
    return rmsd[0]


def gather_pdbfiles (strpath):
    '''

    :param strpath: (str) absolut path to the structure location
    :return: filelist: (list) a list of str

    '''

    filelist = []
    for myfile in os.listdir(strpath):
        if myfile.endswith('.pdb'):
            filelist.append(myfile)
    return filelist

if __name__ == '__main__':
    #pymol.pymol_argv = ['pymol', '-qc'] + sys.argv[1:]
    #pymol.finish_launching()
    #cmd = pymol.cmd

    #test1 = '65-79+87-93+110+111+150+169-180+213-235+290-293+368-371+395-401+420-4220'
    #TODO take the selection string from a parser instead than hardcoded
    test2 = '69-75+368'
    #if len(sys.argv) != 3:
    #    print("Usage:\n gen_rmsd_matrix.py path_to_structures align_string\n")
    #    sys.exit()


    #mobpdb = sys.argv[1]
    #refpdb = sys.argv[2]

    #DEBUGING
    strpath = sys.argv[1]
    #align_str = sys.argv[2]
    align_str = test2

    pdbfiles = []
    if os.path.isabs(strpath) == False:
        try:
            temp = os.path.abspath(strpath)
            strpath = temp
            pdbfiles = gather_pdbfiles(strpath)
        except OSError:
            message = "ERROR: Path %s not found" % strpath
            print(message)
            message = "Please provide a valid path to where the pdb files to process are located"
            print(message)
            sys.exit()
    else:
        try:
            pdbfiles = gather_pdbfiles(strpath)
        except:
            message = "ERROR: Path %s not found" % strpath
            print(message)
            message = "Please provide a valid path to where the pdb files to process are located"
            print(message)
            sys.exit()

    #print pdbfiles
    os.chdir(strpath)

    rmsd_matrix = {}
    i = 0

    for ref in pdbfiles:
        print("=========================================")
        print("Processing pdb: %s as reference" % ref)
        j = 0
        rmsd_matrix[i] = {}
        for mobile in pdbfiles:
            if mobile == ref:
                rmsd_matrix[i][j] = (ref, mobile, "0.00")
            else:
                print("Mobile: " + mobile)
                rmsd = calc_rmsd(mobile, ref, align_str)
                two_rmsd = "%.2f" % rmsd

                if rmsd > 0 and rmsd < 30:
                    rmsd_matrix[i][j] = (ref, mobile, two_rmsd)
                else:
                    print("WARNING: Unresonable RMSD value of %f found, the resulting matrix will likely fail during spectral clustering" % rmsd)
                    print("WARNING: A value of 30 A has been assigned to attempt to savage the RMSD matrix. Be cautious")
                    print("TIP: Exclude structures %s and %s and repeat the calculation" % (ref, mobile))
                    rmsd_matrix[i][j] = (ref, mobile, 30.0)
            j = j + 1
        i = i + 1

    with open('rmsd_matrix.pickle', 'w') as f:
        pickle.dump(rmsd_matrix, f)

    #calc_rmsd(pdbfiles[0],pdbfiles[1], test2)
    #calc_rmsd('2p4j_A.pdb','1sgz_A.pdb', test1)
    #calc_rmsd('2vj7_A.pdb', '3dv1_A.pdb', test1)
    #print("##################################")
    #calc_rmsd('2p4j_A.pdb', '1sgz_A.pdb', test2)
    #calc_rmsd('2vj7_A.pdb', '3dv1_A.pdb', test2)

