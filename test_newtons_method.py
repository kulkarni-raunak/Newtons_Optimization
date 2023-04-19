import pytest
from newtons_method import optimisation, F
import numpy as np

# Test the optimisation class using the given function F
def test_optimisation_class():
    # Define the initial values
    x_initial = np.array([[10], [100], [0]])
    threshold = 10e-6
    num_of_iterations = 100
    # Initialize the optimization class
    problem_00 = optimisation(x_initial, threshold, num_of_iterations, F)
    # Check if the optimization result is close enough to the expected result
    expected_result = np.array([[1.75717907], [4.66261117], [3.82115468]])
    result, n_iter = problem_00.newton_method()
    assert np.allclose(result, expected_result, rtol=1e-3, atol=1e-3)
