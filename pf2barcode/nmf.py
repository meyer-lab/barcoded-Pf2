import numpy as np
from . import matrix_utils as mu
from tqdm import tqdm
from .nnls import nnlsm_blockpivot
from sklearn.decomposition._nmf import _initialize_nmf  # type: ignore


class NMF_ANLS_BLOCKPIVOT:
    """NMF algorithm: ANLS with block principal pivoting

    J. Kim and H. Park, Fast nonnegative matrix factorization: An active-set-like method and comparisons,
    SIAM Journal on Scientific Computing,
    vol. 33, no. 6, pp. 3261-3281, 2011.
    """

    default_max_iter = 100
    default_max_time = np.inf

    def __init__(self, default_max_iter=50, default_max_time=np.inf):
        self.set_default(default_max_iter, default_max_time)

    def set_default(self, default_max_iter, default_max_time):
        self.default_max_iter = default_max_iter
        self.default_max_time = default_max_time

    def run(self, A, k, init=None, max_iter=None):
        """Run a NMF algorithm

        Parameters
        ----------
        A : numpy.array or scipy.sparse matrix, shape (m,n)
        k : int - target lower rank

        Optional Parameters
        -------------------
        init : (W_init, H_init) where
                    W_init is numpy.array of shape (m,k) and
                    H_init is numpy.array of shape (n,k).
                    If provided, these values are used as initial values for NMF iterations.
        max_iter : int - maximum number of iterations.
                    If not provided, default maximum for each algorithm is used.

        Returns
        -------
        (W, H, rec)
        W : Obtained factor matrix, shape (m,k)
        H : Obtained coefficient matrix, shape (n,k)
        rec : dict - (debugging/experimental purpose) Auxiliary information about the execution
        """
        norm_A = mu.norm_fro(A)

        W, H = _initialize_nmf(A, k, init=init)
        H = H.T

        tq = tqdm(range(1, max_iter if max_iter is not None else self.default_max_iter))
        for i in tq:
            # algorithm-specific iteration solver
            (W, H) = self.iter_solver(A, W, H, k, i)

            rel_error = mu.norm_fro_err(A, W, H, norm_A) / norm_A
            tq.set_postfix(rel_error=rel_error, refresh=False)

        W, H, weights = mu.normalize_column_pair(W, H)

        rec = None

        return (W, H, rec)

    def iter_solver(self, A, W, H, k, it):
        Sol, info = nnlsm_blockpivot(W, A, init=H.T)
        H = Sol.T
        Sol, info = nnlsm_blockpivot(H, A.T, init=W.T)
        W = Sol.T
        return (W, H)
