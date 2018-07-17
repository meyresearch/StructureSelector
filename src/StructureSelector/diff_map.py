"""
Implements Diffusion maps on RMSD matrix of pdbs
Author: Antonia Mey <antonia.mey@gmail.com>
"""
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsl

class DiffMap(object):
    """
    Simple implementation of a diffusion map
    """
    def __init__(self, distance_matrix):
        self.K = distance_matrix
        self.diffusion_kernel = None
        self.epsilon = None

    def compute_diffusion_kernel(self, epsilon=None, kernel='gaussian'):
        r""" Computes the Kerneal of a distance matrix

        Parameters:
        -----------

        Returns:
        -------
        """
        if self.epsilon is None:
            self.epsilon = epsilon
        if self.epsilon is None:
            raise("Error: no epsilon given")

        self.diffusion_kernel = np.exp(-self.K**2/self.epsilon)
        else:
            raise("Error: Kernel type not understood.")
        
    def compute_eval_evec_diffusion_process(self):
        r""" Computes the n dominanat eigenvalues and eigenvectors of the diffusion process

        Parameters:
        -----------

        Returns:
        -------
        """

        q, right_norm_vec = self._get_right_norm_vec(weights)
        P = self.normalise_kernel(right_norm_vec)
        evecs, evals = self.get_diffusion_info(P)

    def normalise_kernel(self, right_norm_vec):
        # Perform right normalization
        m = right_norm_vec.shape[0]
        Dalpha = sps.spdiags(right_norm_vec, 0, m, m)
        kernel_matrix = self.kernel_matrix * Dalpha

        # Perform  row (or left) normalization
        row_sum = kernel_matrix.sum(axis=1).transpose()
        n = row_sum.shape[1]
        Dalpha = sps.spdiags(np.power(row_sum, -1), 0, n, n)
        P = Dalpha * kernel_matrix
        return P

    def get_diffusion_info(self):
        raise("Error: not implemented yet")

    def _get_right_norm_vec(self, weights):
        raise("Error: not implemented yet")
