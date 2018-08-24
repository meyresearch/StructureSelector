
import sys, os
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist




def pdbtofasta(pdbfile, three2one):
    '''

    :param pdbfile: (file opened with reading permissions)
    :param three2one: (dictionary containing the three letter to one letter aminoacid code
    :return: str (FASTA sequnce like str)

    Transforms the PDB sequence into a FASTA like str

    '''

    fastaseq = ''
    while True:
        line = pdbfile.readline()
        mor = line.split()
        if line == '':
            break
        elif mor[0] == "ATOM" and len(mor) > 4 and mor[2] == "CA":
            aminoacid = mor[3]
            if len(aminoacid) > 3:
                newamino = aminoacid.lstrip('A')
                if len(newamino) ==3:
                    print("Non standard 4 letter code aminoacid %s found. Assuming multiple occupancy label." % aminoacid)
                    aminoacid = newamino
                else:
                    print("Ignoring residue with 4 letter code %s. Assuming multiple occupancy label " % aminoacid)
            try:
                fastaseq = fastaseq + three2one[aminoacid]
            except KeyError:
                print("Residue %s is not a standard residue and will not be added to the FASTA sequence " % aminoacid)

    return fastaseq

def readfasta(fastafile, valid_aa):
    '''

    :param fastafile: file opened with reading permissions.
    :param valid_aa: Dictionary of standard residues
    :return: fastaseq: A fasta sequence as a str

    '''
    fastaseq = ""
    try:
        while True:
            line = fastafile.readline().rstrip()
            if line == '':
                break
            elif '>' not in line:
                for aa in line:
                   if aa in valid_aa:
                       fastaseq = fastaseq + aa
                   else:
                       print("Error: Unrecognized character %s. Are you sure this is a valid fasta seq?") % aa
                       raise IOError
    except:
        raise

    return fastaseq

def getalignement(seq1, seq2):
    '''


    :param seq1: A FASTA like str
    :param seq2: A FASTA like str.
    :return: alignment: str

    Uses biopython pairwise function to align both sequences. Since it is intended to match identical sequences the default
    score parameters:

    Match score: +5
    Mismatch score: -5
    NewGaps: -10
    Extend a Gap = -1


    '''
    #matrix = matlist.blosum62
    alignment = pairwise2.align.globalms(seq1, seq2, 5,-5,-10,-1)
    #print alignment[0][0]
    return alignment[0][0]

def rewrite_pdb_with_FASTA_num(pdbfile, alignment):
    '''

    :param pdbfile:
    :param alignment:
    :return:

    '''

    new_index = calculate_first_residue(alignment)
    #new_residue = new_index + 1 # This correct the 0 to 1 first index between python and PDB
    new_residue = new_index
    # These keep track of the residues in the PDB
    last_residue = 0
    current_residue = -99
    line_tracker = 0

    with open('_tempfile', 'w') as outf:
        while True:
            line = pdbfile.readline()
            line_tracker = line_tracker + 1
            if line == '':
                break
            elif "ATOM" in line:
                mor=split_ATOM_PDB_line(line)
                # Check if it is a new residue
                if int(mor[6]) !=  current_residue:

                    # update the residue control vars
                    last_residue = current_residue
                    current_residue = int(mor[6])

                    residue_checkpoint = new_index # track the last residue (for missing C-terminal regions)

                    while True:
                        # Ensure new_index is not bigger thant the sequence lenght (for missing  C-terminal regions)
                        if new_index > len(alignment) - 1 :
                            new_residue = residue_checkpoint + 1
                            break
                        elif alignment[new_index] == '-':
                            # shift new_residues for as many residues are in the gap.
                            new_residue = new_residue + 1
                            new_index = new_index + 1
                        else:
                            # otherwise increment new_residue and new_index by one
                            new_residue = new_residue + 1
                            new_index = new_index + 1
                            break
                try:
                    newline="{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          " \
                            "{:>2s}{:2s}\n".format(mor[0],int(mor[1]),mor[2],mor[3],mor[4],mor[5],new_residue,mor[7],
                                                 float(mor[8]),float(mor[9]),float(mor[10]),float(mor[11]),float(mor[12]),
                                                 mor[13],mor[14])
                    #outf.write(newline)
                except IndexError:
                    message="Error: line %d does not seem to be in standard pdb format: \n\n %s \n\n" % (line_tracker, line)
                    print(message)
                    message = "Expected a split len of 12, found %d" % len(mor)
                    print(message)
                    raise
                    #sys.exit()
            elif "TER" in line:
                newline = line
            outf.write(newline)

        outf.write("END")

def split_ATOM_PDB_line(pdbline):
    '''

    :param pdbline: a str containing an ATOM entry from a PDB file
    :return: splitline: a list containing a split of an ATOM pdbline folloowing the format outlined below

    Index   Descipton                                       lenght  format      range   string slicing
    1 	    "ATOM " or "HETATM" 	                        6 	    {:6s}   	01-06 	[0:6]
    2 	    atom serial number 	                            5 	    {:5d} 	    07-11 	[6:11]
    3 	    atom name 	                                    4 	    {:^4s} 	    13-16 	[12:16]
    4 	    alternate location indicator 	                1 	    {:1s} 	    17 	    [16:17]
    5 	    residue name 	                                3 	    {:3s} 	    18-20 	[17:20]
    6 	    chain identifier 	                            1 	    {:1s} 	    22 	    [21:22]
    7 	    residue sequence number 	                    4 	    {:4d} 	    23-26 	[22:26]
    8 	    code for insertion of residues 	                1 	    {:1s} 	    27 	    [26:27]
    9 	    orthogonal coordinates for X (in Angstroms) 	8 	    {:8.3f} 	31-38 	[30:38]
    10 	    orthogonal coordinates for Y (in Angstroms) 	8 	    {:8.3f} 	39-46 	[38:46]
    11 	    orthogonal coordinates for Z (in Angstroms) 	8 	    {:8.3f} 	47-54 	[46:54]
    12 	    occupancy 	                                    6 	    {:6.2f} 	55-60 	[54:60]
    13 	    temperature factor 	                            6 	    {:6.2f} 	61-66 	[60:66]
    14 	    element symbol 	                                2 	    {:>2s} 	    77-78 	[76:78]
    15 	    charge on the atom 	                            2 	    {:2s} 	    79-80 	[78:80]

    Source: Pierre Poulain, Universite Paris Diderot - Paris 7
            http://cupnet.net/pdb-format/

    '''

    splitline = []
    splitline.append(pdbline[0:6])
    splitline.append(pdbline[6:11])
    splitline.append(pdbline[12:16])
    splitline.append(pdbline[16:17])
    splitline.append(pdbline[17:20])
    splitline.append(pdbline[21:22])
    splitline.append(pdbline[22:26])
    splitline.append(pdbline[26:27])
    splitline.append(pdbline[30:38])
    splitline.append(pdbline[38:46])
    splitline.append(pdbline[46:54])
    splitline.append(pdbline[54:60])
    splitline.append(pdbline[60:66])
    splitline.append(pdbline[76:78])
    splitline.append(pdbline[78:80])

    return splitline





def calculate_first_residue(alignment):
    '''

    :param alignment: FASTA like sequence where '-' indicate gaps
    :return: firstresidue: int
    '''

    firstresidue = 0
    for res in alignment:
        if res != '-':
            break
        firstresidue = firstresidue + 1

    return firstresidue


if __name__ == '__main__':

    one2threeAA = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'E': 'GLU', 'Q': 'GLN', 'G': 'GLY',
                   'H': 'HIS',
                   'I': 'ILE', 'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR',
                   'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}

    three2oneAA = {}

    for resi in one2threeAA:
        tempval = one2threeAA[resi]
        three2oneAA[tempval] = resi

    if len(sys.argv) != 3:
        message = "Usage: cleanPDB.py PDB FASTA\n PDB: Single Chain PDBfile\n FASTA: Reference aa sequence file"
        print(message)
        sys.exit()




    pdbfile = sys.argv[1]
    fastafile = sys.argv[2]

    message = "Procesing %s..." % pdbfile
    print(message)

    #TODO split between complete and incomplete structures.


    with open(fastafile) as ff:
        ref_fasta = readfasta(ff, one2threeAA)
    message ="The reference FASTA sequence of the Protein is:\n" + ref_fasta
    print(message)
    with open(pdbfile) as pdbf:
        pdbfasta = pdbtofasta(pdbf, three2oneAA)
        pdbf.seek(0)

    alignment = getalignement(pdbfasta, ref_fasta)
    with open(pdbfile) as pdbf:
        rewrite_pdb_with_FASTA_num(pdbf, alignment)
    command = "mv _tempfile %s" % pdbfile
    os.system(command)




    sys.exit()