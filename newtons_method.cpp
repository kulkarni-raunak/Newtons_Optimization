#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>    // std::transform
#include <functional>   // std::plus
#include "helper_function.h"
// #include "helper_function.cpp"
using namespace std;


class Optimisation {        // The class
  public:          // Access specifier
    float threshold;
    int num_iterations;
    float step_size = 10e-6;
    vector<float> initial_value;


    Optimisation(float threshold_, int num_iterations_, vector<float>& initial_value_) {  // Constructor with parameters
      threshold = threshold_;
      num_iterations = num_iterations_;
      initial_value.assign(initial_value_.begin(),initial_value_.end());
    }

    vector<float> func_value(vector<float> initial_value_, vector<float> (*Func)(vector<float>)){
      return Func(initial_value_);
    }

    float dFi_dxj(vector<float> x_value, int i, int j, vector<float> (*Func)(vector<float>)){
      /*
    This function calculates each element of the Jacobian matrix using finite differences. 
    Inputs are 'x' = vector, 'i' = to select i-th function from Function matrix
                and 'j' = indicates j-th element in input vector.
    Output is element of Jacobian matrix.
    */
      vector<float> x_value_stepsize_e_j(sizeof(x_value),0.0);
      x_value_stepsize_e_j.at(j) = step_size;
      std::transform (x_value_stepsize_e_j.begin(), x_value_stepsize_e_j.end(), x_value.begin(), x_value_stepsize_e_j.begin(), std::plus<float>());
      float dFi_dxj_ = (Func(x_value_stepsize_e_j)[i] - Func(x_value)[i])/step_size; // returns (F_i(x + step_size*e_j) - F_i(x))/step_size 
      return dFi_dxj_;
    }

    vector<vector<float>> jacobian_matrix(vector<float> x_value_, vector<float> (*Func)(vector<float>)){
      /*
      This function calculates the Jacobian matrix w.r.t Function matrix and input vector 'x'.
      Input is 'x' = vector.
      Output is a jacobian matrix.
      */
      int i = Func(x_value_).size();
      int j = x_value_.size();
      vector<vector<float>> j_mat(i, vector<float> (j,0));
      for (int i_i = 0; i_i < i; i_i++)
      {
        for (int j_i = 0; j_i < j; j_i++)
        {
          j_mat[i_i][j_i] = dFi_dxj(x_value_, i_i, j_i, (*Func));
        }
      }
      return j_mat;
    }

    std::tuple <vector<float>, int> newton_method(vector<float> (*Func)(vector<float>)){
      /*
      This function performs newton's method of optimisation for a defined function, initial value and specific number of iterations.
      The function returns the minimum value if a difference between estimated values is below a threshold or maximum iterations are reached. 
      */
      int size_i = initial_value.size();
      vector<float> x_old(size_i, 0.0), x_new(size_i, 0.0), delta_x(size_i, 0.0);
      x_old.assign(initial_value.begin(),initial_value.end());
      for(int iter_num=1; iter_num <= num_iterations; iter_num++){

        vector<vector<float>> j_mat =jacobian_matrix(x_old,(*Func));

        vector<vector<float>> inv_j_mat = getInverse(j_mat);


        vector<float> mult_matrix = multiple_mat_vec(inv_j_mat,Func(x_old));

        vector<float> negative_vector = negate_vector(mult_matrix);


        // x_new = x_old - Jacobian_matrix_inverse(x) * F(x)
        std::transform (x_old.begin(), x_old.end(), negative_vector.begin(), x_new.begin(), std::plus<float>());

        vector<float> neg_x_old = negate_vector(x_old);

        // ||x_new-x_old|| 2nd-norm 
        std::transform (x_new.begin(), x_new.end(), neg_x_old.begin(), delta_x.begin(), std::plus<float>());
        float L2_norm = 0.0;
        for (int k = 0; k < delta_x.size(); k++)
        {
          L2_norm += delta_x[k]*delta_x[k];
        }
        L2_norm = sqrt(L2_norm);
        cout<< "L2_norm is " << L2_norm << endl;
        if (L2_norm <= threshold)
        {
          cout << "Reached termination stage" << endl;
          return std::make_tuple(x_new, iter_num);
        }
        x_old.assign(x_new.begin(),x_new.end());
      }
      cout << "Reached last iteration" << endl;
      return std::make_tuple(x_old, num_iterations);      
    }

};

vector<float> Func(vector<float> x){
  /*
  Function defined for input vector x = x1,x2,x3
  F1 = x1/x2 + x3/x1
  F2 = (0.5 * x2 ^ 3) -(250 * x2 * x3) - (75000 * x3 ^ 2)
  F3 = (e ^ -x3) + (x3 * e)

  Input is  3-D vector 'x'
  Output is 3-D vector
  */
  vector<float> F;
  F.push_back((x[0]/x[1]) + (x[2]/x[0]));
  F.push_back((0.5 * pow(x[1],3)) - (250 * x[1] *x[2]) - (75000 * pow(x[2], 2)));
  F.push_back(exp(-1 * x[2]) + (x[2] * exp(1)));
  return F;
}

int main() {

  // Problem setup for optimization
  float x_initial[3] = {10, 100, 0};
  vector<float> initial_value(x_initial, x_initial + sizeof(x_initial) / sizeof(float) );

  Optimisation Problem_00(10e-6, 100, initial_value);
  vector<float> func_value = Problem_00.func_value(initial_value, &Func);
  cout << "Completed problem initialisation" << endl;

  // Finding minimum using newtons method
  vector<float> final_value;
  int iterations_completed;
  tie(final_value, iterations_completed) = Problem_00.newton_method(&Func);


  // Print values
  cout << "F_minimum is found at " << final_value.at(0) << " " << final_value.at(1) << " " << final_value.at(2) << endl;
  cout << "Completed in " << iterations_completed << endl;


  return 0;
}
