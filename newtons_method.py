# import necessary libraries 
import numpy as np 
        
class optimisation:
  """
  Optimisation class currently only capable to optimise using newtons method.
  """

  def __init__(self, initial_value, threshold, num_ierations, Function_):
    """
    Initialise all the attributes
    """
    self.initial_value = initial_value
    self.threshold = threshold
    self.num_ierations = num_ierations
    self.Function_ = Function_
    self.step_size = 10e-6
  
  def dFi_dxj(self, x, i, j):
    """
    This function calculates each element of the Jacobian matrix using finite differences. 
    Inputs are 'x' = vector, 'i' = to select i-th function from Function matrix
                and 'j' = indicates j-th element in input vector.
    Output is element of Jacobian matrix. 
    """
    e_j = np.zeros((len(x),1)) # Null vector same in dimensions as input vector.
    e_j[j] = 1 # Turns null vector in unit vector at j-th position.
    dFi_dxj_ = (self.Function_(x + self.step_size * e_j)[i]
                - self.Function_(x)[i])/self.step_size # returns (F_i(x + step_size*e_j) - F_i(x))/step_size 
    return dFi_dxj_

  def jacobian_matrix(self, x):
    """
    This function calculates the Jacobian matrix w.r.t Function matrix and input vector 'x'.
    Input is 'x' = vector.
    Output is a jacobian matrix.
    """
    i, j = len(self.Function_(x)),len(x)
    j_mat = [[self.dFi_dxj(x, i_i, j_i)[0] for j_i in range(j)] for i_i in range(i)]
    return j_mat

  def newton_method(self):
    """
    This function performs newton's method of optimisation for a defined function, initial value and specific number of iterations.
    The function returns the minimum value if a difference between estimated values is below a threshold or maximum iterations are reached. 
    """
    x_old = self.initial_value
    for iter_num in range(self.num_ierations):
        # x_new = x_old - Jacobian_matrix_inverse(x) * F(x)
        x_new = x_old - np.matmul(np.linalg.inv(self.jacobian_matrix(x_old)), self.Function_(x_old)) 
        delta_x = np.linalg.norm(x_old - x_new, 2) # ||x_new-x_old|| 2nd-norm 
        if delta_x <= self.threshold:
            return x_new, (iter_num+1)
        x_old = x_new
    return x_old, (iter_num+1)


def F(x):
    """
    Function defined for input vector x = x1,x2,x3
    F1 = x1/x2 + x3/x1
    F2 = (0.5 * x2 ^ 3) -(250 * x2 * x3) - (75000 * x3 ^ 2)
    F3 = (e ^ -x3) + (x3 * e)

    Input is  3-D vector 'x'
    Output is 3-D vector
    """
    assert len(x) == 3 # Check if the input vector is compactible with the function

    F1 = (x[0]/x[1]) + (x[2]/x[0])
    F2 = (0.5 * x[1] ** 3) - (250 * x[1] * x[2]) - (75000 * x[2] ** 2)
    F3 = np.exp(-1 * x[2]) + (x[2] * np.exp(1))
    return np.array([F1, F2, F3])


#  Problem setup for optimization
x_initial = np.array([[10], [100], [0]])
threshold = 10e-6
num_of_iterations = 100
problem_00 = optimisation(x_initial, threshold, num_of_iterations, F)

# Finding minimum using newtons method
x_prime, n_iter = problem_00.newton_method()

print(f"X_min is {x_prime} and {n_iter} iterations completed")