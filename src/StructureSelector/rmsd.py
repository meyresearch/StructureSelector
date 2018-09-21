import pymol
import pickle
import sys


class RMSD(object):
    def __init__(self, mobile_pdb_list, align_string='69-75+368', reference_pdb=None):
        self.matrix = {}
        self.reference_pdb = reference_pdb
        self.mobile_pdb_list = mobile_pdb_list
        self.align_str = align_string

    def compute_rmsd_matrix(self):
        """
        computes the RMSD matrix
        """

        i = 0
        for ref in self.mobile_pdb_list:
            print("=========================================")
            print("Processing pdb: %s as reference" % ref)
            j = 0
            self.matrix[i] = {}
            for mobile in self.mobile_pdb_list:
                if mobile == ref:
                    self.matrix[i][j] = (ref, mobile, "0.00")
                else:
                    print("Mobile: " + mobile)
                    rmsd = self._calc_rmsd(mobile, ref, self.align_str)
                    two_rmsd = "%.2f" % rmsd

                    if rmsd > 0 and rmsd < 30:
                        self.matrix[i][j] = (ref, mobile, two_rmsd)
                    else:
                        print(
                            "WARNING: Unresonable RMSD value of %f found, the resulting matrix will likely fail during spectral clustering" % rmsd)
                        print(
                            "WARNING: A value of 30 A has been assigned to attempt to savage the RMSD matrix. Be cautious")
                        print("TIP: Exclude structures %s and %s and repeat the calculation" % (ref, mobile))
                        self.matrix[i][j] = (ref, mobile, 30.0)
                j = j + 1
            i = i + 1

    def write_matrix(self, filename='rmsd_matrix.pickl'):
        with open(filename, 'w') as f:
            pickle.dump(self.matrix, f)

    def _calc_rmsd(self, mobpdb, refpdb, align_str):
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
                         "resi %s and (%s)" % (align_str, ref_obj), transform=0, object="test4", reset=1)

        cmd.delete(mobile_obj)
        cmd.delete(ref_obj)

        # print(ca_rmsd)
        # print(rmsd)
        return rmsd[0]

    #def align_pdb(self, mobile_pdb, reference_pdb):