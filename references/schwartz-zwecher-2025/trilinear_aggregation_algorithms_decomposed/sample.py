import numpy as np


def main():
    # load decomposed algorithm from the .npz file
    algorithm = np.load("algorithm_44_44_44_36110_decomposed.npz")
    u_phi, v_phi, w_phi, phi = algorithm["u_phi"], algorithm["v_phi"], algorithm["w_phi"], algorithm["phi"]
    u, v, w = u_phi @ phi, v_phi @ phi, w_phi @ phi

    # generate random matrices
    a = np.random.rand(44, 44)
    b = np.random.rand(44, 44)

    # calculate expected result
    expected_result = a @ b

    # calculate actual result
    actual_result = w.T @ ((u @ a.ravel()) * (v @ b.ravel()))
    actual_result = actual_result.reshape((44, 44)).T

    # check if the results are equal (within tolerance)
    if np.allclose(expected_result, actual_result):
        print("The results match.")
    else:
        print("The results do not match.")


if __name__ == '__main__':
    main()
