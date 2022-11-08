import numpy as np


def svd_compress(kspace, out_coils=None, tol=0.1):
    """
    Simple implementation of singular value decomposition (SVD) coil compression
    Parameters
    ---------
    kspace : numpy array
    out_coils: int
        Number of output coils
    tol:
        Threshold (if number of output coils is not given)
    Returns
    -------
    numpy array
        Coil compressed k-space

    If out_coils is none,
    """
    X = kspace.reshape(-1, kspace.shape[-1])
    _, S, Vh = np.linalg.svd(X, full_matrices=False)
    V = Vh.T.conj()
    if not out_coils:
        out_coils = sum(s > tol for s in S / S.max())
    compressed = X @ V[:, :out_coils]
    return compressed.reshape(kspace.shape[:-1] + (out_coils,))

# https://stackoverflow.com/questions/19636487/concatenation-of-numpy-arrays-of-unknown-dimension-along-arbitrary-axis/19639323
def atleast_nd(x, n):
    return np.array(x, ndmin=n, subok=True, copy=False)

# From MRSRL (https://github.com/MRSRL/mridata-recon/blob/master/recon_2d_fse.py)
def isrmrmd_user_param_to_dict(header):
    """
    Store ISMRMRD header user parameters in a dictionary.
    Parameter
    ---------
    header : ismrmrd.xsd.ismrmrdHeader
        ISMRMRD header object
    Returns
    -------
    dict
        Dictionary containing custom user parameters
    """
    user_long = list(header.userParameters.userParameterLong)
    user_double = list(header.userParameters.userParameterDouble)
    user_string = list(header.userParameters.userParameterString)
    user_base64 = list(header.userParameters.userParameterBase64)

    return {
        entry.name: entry.value
        for entry in user_long + user_double + user_string + user_base64
    }

# From sigpy (https://github.com/mikgroup/sigpy/blob/master/sigpy/util.py)
def resize(ary, oshape, ishift=None, oshift=None):
    """Resize with zero-padding or cropping.
    Args:
        ary (array): Input array.
        oshape (tuple of ints): Output shape.
        ishift (None or tuple of ints): Input shift.
        oshift (None or tuple of ints): Output shift.
    Returns:
        array: Zero-padded or cropped result.
    """

    ishape1, oshape1 = _expand_shapes(ary.shape, oshape)

    if ishape1 == oshape1:
        return ary.reshape(oshape)

    if ishift is None:
        ishift = [max(i // 2 - o // 2, 0) for i, o in zip(ishape1, oshape1)]

    if oshift is None:
        oshift = [max(o // 2 - i // 2, 0) for i, o in zip(ishape1, oshape1)]

    copy_shape = [min(i - si, o - so) for i, si, o, so in zip(ishape1, ishift, oshape1, oshift)]
    islice = tuple(slice(si, si + c) for si, c in zip(ishift, copy_shape))
    oslice = tuple(slice(so, so + c) for so, c in zip(oshift, copy_shape))

    output = np.zeros(oshape1, dtype=ary.dtype)
    ary = ary.reshape(ishape1)
    output[oslice] = ary[islice]

    return output.reshape(oshape)


# From sigpy (https://github.com/mikgroup/sigpy/blob/master/sigpy/util.py)
def _expand_shapes(*shapes):

    shapes = [list(shape) for shape in shapes]
    max_ndim = max(len(shape) for shape in shapes)
    shapes_exp = [[1] * (max_ndim - len(shape)) + shape for shape in shapes]

    return tuple(shapes_exp)