import pymol
import pickle
import os


class RMSD(object):
    def __init__(self, mobile_pdb_list, align_string=None, verbose = False):
        """

        Parameters:
        ----------
        mobile_pdb_list : list
            list containing path to all pdbfiles used for the RMSD matrix
        align_string : string
            pymol string used for alignment e.g. 140-150+200, means alignment should be done on residues 140 to 150 and 200
        verbose: boolean
            Be loud and noisy
        """
        self.matrix = {}
        self.mobile_pdb_list = mobile_pdb_list
        self.align_str = align_string
        self.rmsd_align = False
        if align_string == None:

            self.rmsd_align = True
        self.verbose = verbose


    def compute_rmsd_matrix(self):
        """
        computes the RMSD matrix between all pdbs
        TODO: add checks that all the proteins are indeed the same and are alignable
        """
        if self.rmsd_align:
            print('No selection string was passed therefore, RMSD distances and alignment will be carried out'
                  'on all C-alpha atoms of the protein')
        else:
            print('Processing RMSD computation for %d proteins using the selection string: resid %s' %(len(self.mobile_pdb_list), self.align_str))

        i = 0
        for ref in self.mobile_pdb_list:
            if self.verbose:
                print("=========================================")
                print("Processing pdb: %s as reference" % ref)
            ref_base = os.path.basename(ref)
            j = 0
            self.matrix[i] = {}
            for mobile in self.mobile_pdb_list:
                mobile_base = os.path.basename((mobile))
                if mobile == ref:
                    self.matrix[i][j] = (ref_base, mobile_base, "0.00")
                else:
                    if self.verbose:
                        print("Mobile: " + mobile)
                    rmsd = self._calc_rmsd(mobile, ref)
                    two_rmsd = "%.2f" % rmsd

                    if rmsd > 0 and rmsd < 30:
                        self.matrix[i][j] = (ref_base, mobile_base, two_rmsd)
                    else:
                        print(
                            "WARNING: Unresonable RMSD value of %f found, the resulting matrix will likely fail during spectral clustering" % rmsd)
                        print(
                            "WARNING: A value of 30 A has been assigned to attempt to savage the RMSD matrix. Be cautious")
                        print("TIP: Exclude structures %s and %s and repeat the calculation" % (ref, mobile))
                        self.matrix[i][j] = (ref_base, mobile_base, 30.0)
                j = j + 1
            i = i + 1

    def write_matrix(self, filename='rmsd_matrix.pickl'):
        """Writes out a pickel file with the RMSD matrix

        Parameters:
        -----------
            filename: string
                Filename of the pickel object that contains the RMSD matrix dictionary
        """
        with open(filename, 'wb') as f:
            pickle.dump(self.matrix, f)

    def _calc_rmsd(self, mobpdb, refpdb):
        """ Uses pymol functions to calculate the rmsd among atoms in align_str between the mobpdb and refpdb

        Parameters:
        -----------
        modpdb : string
            string with path to the mobile pdb file
        refpdb : string
            string with path to the reference pdb file

        Returns:
        --------
        rmsd : float
            vlaue of RMSD of the two pdb files.
        """

        pymol.pymol_argv = ['pymol', '-qc']
        pymol.finish_launching()
        cmd = pymol.cmd

        mobile_obj = 'mobile'
        ref_obj = 'reference'

        cmd.load(mobpdb, mobile_obj)
        cmd.load(refpdb, ref_obj)


        rmsd = None
        if self.rmsd_align:
            rmsd = cmd.align("polymer and name CA and (%s)" % mobile_obj, "polymer and name CA and (%s)" % ref_obj,
                            object="test4", reset=0)
        else:
            rmsd = cmd.align("resi %s and (%s)" % (self.align_str, mobile_obj),
                         "resi %s and (%s)" % (self.align_str, ref_obj), transform=0, object="test4", reset=1)

        cmd.delete(mobile_obj)
        cmd.delete(ref_obj)

        # print(ca_rmsd)
        if self.verbose:
            print(rmsd)
        return rmsd[0]

    #TODO: Write function that allows to write out align pdbs.
    #def align_pdb(self, mobile_pdb, reference_pdb):