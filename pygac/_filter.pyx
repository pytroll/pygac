import numpy as np
cimport numpy as np

DTYPE = np.double
ctypedef np.double_t DTYPE_T

def _mean_filter(np.ndarray[DTYPE_T, ndim=2] data, int box_size, 
                 DTYPE_T fill_value):
    """Filter a 2D array using an arithmetic mean kernel.

    Compute the arithmetic mean of the valid elements within a box of size
    (boxsize x boxsize) around each pixel. Fill values are not taken into
    account.

    Args:
        data (np.ndarray): 2D array to be filtered. Masked arrays are not
            supported, invalid data must be filled with fill_value.
        box_size (int): Specifies the box_size. Must be odd.
        fill_value: Value indicating missing data.

    Returns:
        np.ndarray: The filtered array
    """
    assert data.dtype == DTYPE

    cdef int nrows = data.shape[0]
    cdef int ncols = data.shape[1]
    cdef int row, col, brow, bcol, num_valid
    cdef DTYPE_T sum_valid
    cdef np.ndarray[DTYPE_T, ndim=2] filtered = fill_value*np.ones(
        (nrows, ncols), dtype=DTYPE)
    cdef int radius = box_size//2

    for row in range(nrows):
        for col in range(ncols):
            # Reset values
            sum_valid = 0.
            num_valid = 0

            # Compute box mean
            for brow in range(max(row-radius, 0), 
                              min(row+radius+1, nrows)):
                for bcol in range(max(col-radius, 0), 
                                  min(col+radius+1, ncols)):
                    if data[brow, bcol] != fill_value:
                        sum_valid += data[brow, bcol]
                        num_valid += 1
            if num_valid > 0:
                filtered[row, col] = sum_valid/float(num_valid)

    return filtered 
