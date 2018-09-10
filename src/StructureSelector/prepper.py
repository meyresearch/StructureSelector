"""
Preps a list of pdb files for computing RMSD matrix
Author: Antonia Mey <antonia.mey@gmail.com>
Author: Jordi Juarez Jimenez
"""
import requests
import os

class Prepper(object):

    errorcode = 99

    def download_pdbs(self, pdb_list, directory='/home/ppxasjsm/Downloads/', link='http://files.rcsb.org/download/'):
        r""" Automatically downloads pdbs from a list of 4 character pdb identifiers

        Parameters:
        -----------
        pdb_list : list
            list containing names of pdb strings
        directory : string
            string to the download directory
        link : string

        Returns:
        success : int
            error code indicating a successful or failed download
        """

        for pdb in pdb_list:
            pdb_ID = pdb.strip().upper()
            pdb_file = pdb_ID+'.pdb'
            download_file = os.path.join(directory, pdb_file)
            print(pdb_file)
            if not os.path.isfile(download_file):
                print('File does not exist retrieving pdbid %s' %pdb_ID)
                url = link+pdb_file
                r = requests.get(url)
                if r:
                    f = open(download_file, 'bw')
                    f.write(r.content)
                    f.close()


    def clean_pdb(self, pdbfile):

        """ Cleans input pdb file by deleting all entries differen than ATOM TER and END
        Parameters:
        ----------
        pdbfile : filehandle
            name of the pdbfile to be worked on



        Returns:
        -------
        errorcode : int
            success or failure in cleaning pdb
        """

        try:
            with open('tempfile.pdb', 'w') as f:

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
            return self.errorcode

    def split_chains(self, pdbfile, outrootname):

        """ Rewrites pdbfile splitting it in one file for each chain
        Parameters:
        -----------
        pdbfile : filehandle
            pdbfile that contains multiple chains to be split

        Returns:
        -------
        errorcode : int
            success of failure in splitting chains
        """

        allchains = {}
        currentchain = ""
        try:
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
                with open(outrootname + "_%s.pdb" % chain, 'w') as f:
                    for line in allchains[chain]:
                        f.write(line)

            return 0

        except IOError:
            raise
            return self.errorcode

    def pdb_to_fasta(self, pdbfile, three2one):

        """Transforms the PDB sequence into a FASTA like str
        Parameters:
        -----------
        pdbfile : file
            pdbfile used for the conversion
        thee2one : dictionary
            mapping 3 letter amino acid codes to one letter codes
        Returns:
        --------
        fastaseq : str
            FASTA sequence like string
        """

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
                    if len(newamino) == 3:
                        print(
                            "Non standard 4 letter code aminoacid %s found. Assuming multiple occupancy label." % aminoacid)
                        aminoacid = newamino
                    else:
                        print(
                            "Ignoring residue with 4 letter code %s. Assuming multiple occupancy label " % aminoacid)
                try:
                    fastaseq = fastaseq + three2one[aminoacid]
                except KeyError:
                    print(
                        "Residue %s is not a standard residue and will not be added to the FASTA sequence " % aminoacid)

        return fastaseq
