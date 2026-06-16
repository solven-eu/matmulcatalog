We provide the encoding U, V and decoding W matrices for the algorithms described in the paper for base cases 20
through 50, as defined in our paper "Towards Faster Feasible Matrix Multiplication by Trilinear Aggregation".
The matrices are provided as .npz files, which can be loaded using numpy.
The files are named "algorithm_{n0}_{n0}_{n0}_{number_of_multiplications}_decomposed.npz"
We further provide sample code for testing the algorithms on random inputs in python.

The U, V, W encoding and decoding matrices may be obtained from the decomposed algorithm by computing
U = u_phi * phi
V = v_phi * phi
W = w_phi * phi

The columns order of U, V, W assumes that the inputs matrices A and B and the output matrix C are represented by vectors
in a row-major order. For example, an input
A = [[1, 2],
     [3, 4]]
is represented as the vector [1, 2, 3, 4]^T.

The output matrix is (AB)^T rather than AB.

The correctness of the algorithms can be verified either by generating random inputs and verifying that the output of
the algorithm is correct, or deterministically by using the triple product condition of Brent (1970).
This requires that for all i1, i2, j1, j2, k1, k2, the following equation holds:

<U[:, k2 * n + i1], V[:, i2 * n + j1], W[:, j2 * n + k1]> = delta_{i1,i2} * delta_{j1,j2} * delta_{k1,k2}

where U, V, W are the encoding matrices, and delta_{i,j} is the Kronecker delta function.


