import StructureSelector.prepper as p
import StructureSelector.rmsd as rmsd
import glob
import os
import sys

catS_prep = p.Prepper()

f = open('../../tests/Cathepsin/CatS.csv', 'r')
lines = f.readlines()
f.close()
pdb_list = lines[0].strip().split(', ')

#pdb_files = glob.glob('../../tests/Cathepsin/PDBS/*')
"""if pdb_files is None:
    catS_prep.download_pdbs(pdb_list, directory='../../tests/Cathepsin/PDBS')
    pdb_files = glob.glob('../../tests/Cathepsin/PDBS/*')

for filename in pdb_files:
    outrootname, extension = os.path.splitext(filename)
    print(filename)
    print(outrootname)
    with open(filename) as inpf:
        status = catS_prep.clean_pdb(inpf)
        if status != 0:
            message = "Unfortunately cleanPDB failed to process the following pdb file: %s" % filename
            print(message)
            sys.exit(status)
        else:
            with open('tempfile.pdb') as tmpfile:
                status = catS_prep.split_chains(tmpfile, outrootname)

            command = "rm _tempfile.pdb %s" % filename
            os.system(command)
            message = "Done!"
            print(message)
            #sys.exit(status)
"""
"""fasta_file = '../../tests/Cathepsin/CathS.fasta'
with open(fasta_file) as ff:
    ref_fasta = catS_prep.read_fasta(ff)
    message ="The reference FASTA sequence of the Protein is:\n" + ref_fasta
    print(message)

pdbfile = '../../tests/Cathepsin/2G6D_A.pdb'
with open(pdbfile) as pdbf:
    pdbfasta = catS_prep.pdb_to_fasta(pdbf)
    pdbf.seek(0)

alignment = catS_prep.get_alignement(pdbfasta, ref_fasta)
with open(pdbfile) as pdbf:
    catS_prep.rewrite_pdb_with_FASTA_num(pdbf, alignment)
command = "mv _tempfile %s" % pdbfile
"""

chainA = glob.glob('../../tests/Cathepsin/PDBS/chainA/*.pdb')

r = rmsd.RMSD(chainA)
r.compute_rmsd_matrix()

print(r.matrix)
