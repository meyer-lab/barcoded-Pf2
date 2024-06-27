import numpy as np
import scipy.sparse as sps
import numpy.linalg as nla
import math


def norm_fro(X):
    """Compute the Frobenius norm of a matrix

    Parameters
    ----------
    X : numpy.array or scipy.sparse matrix

    Returns
    -------
    float
    """
    if sps.issparse(X):  # scipy.sparse array
        return math.sqrt(X.multiply(X).sum())
    else:  # numpy array
        return nla.norm(X)


def norm_fro_err(X, W, H, norm_X):
    """Compute the approximation error in Frobeinus norm

    norm(X - W.dot(H.T)) is efficiently computed based on trace() expansion
    when W and H are thin.

    Parameters
    ----------
    X : numpy.array or scipy.sparse matrix, shape (m,n)
    W : numpy.array, shape (m,k)
    H : numpy.array, shape (n,k)
    norm_X : precomputed norm of X

    Returns
    -------
    float
    """
    sum_squared = (
        norm_X * norm_X
        - 2 * np.trace(H.T.dot(X.T.dot(W)))
        + np.trace((W.T.dot(W)).dot(H.T.dot(H)))
    )
    return math.sqrt(np.maximum(sum_squared, 0))


def column_norm(X):
    """Compute the norms of each column of a given matrix

    Parameters
    ----------
    X : numpy.array or scipy.sparse matrix

    Returns
    -------
    numpy.array
    """
    if sps.issparse(X):
        norm_vec = np.sqrt(X.multiply(X).sum(axis=0))
        return np.asarray(norm_vec)[0]
    else:
        norm_vec = np.sqrt(np.sum(X * X, axis=0))
        return norm_vec


def normalize_column_pair(W, H):
    """Column normalization for a matrix pair

    Scale the columns of W and H so that the columns of W have unit norms and
    the product W.dot(H.T) remains the same.  The normalizing coefficients are
    also returned.

    Side Effect
    -----------
    W and H given as input are changed and returned.

    Parameters
    ----------
    W : numpy.array, shape (m,k)
    H : numpy.array, shape (n,k)

    Returns
    -------
    ( W, H, weights )
    W, H : normalized matrix pair
    weights : numpy.array, shape k
    """
    norms = column_norm(W)

    toNormalize = norms > 0
    W[:, toNormalize] = W[:, toNormalize] / norms[toNormalize]
    H[:, toNormalize] = H[:, toNormalize] * norms[toNormalize]

    weights = np.ones(norms.shape)
    weights[toNormalize] = norms[toNormalize]
    return (W, H, weights)
