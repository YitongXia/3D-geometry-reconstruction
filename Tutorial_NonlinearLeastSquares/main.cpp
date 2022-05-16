/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <easy3d/optimizer/optimizer_lm.h>

using namespace easy3d;



/// To use the Levenberg-Marquardt method to solve a non-linear least squares method, we need to define our own
/// objective function that inherits 'Objective_LM'.

class MyObjective : public Objective_LM {
public:
    MyObjective(int num_func, int num_var) : Objective_LM(num_func, num_var) {}

    /**
     *  Calculate the values of each function at x and return the function values as a vector in fvec.
     *  @param  x           The current values of variables.
     *  @param  fvec        Return the value vector of all the functions.
     *  @return Return a negative value to terminate.
     *
     *  NOTE: This function implements f = (x0 - 1.0)^2 + (x1 - 1.0)^2. A client problem must implement
     *      this function to evaluate the values of each function in the expression of x.
     */
    int evaluate(const double *x, double *fvec) {
        fvec[0] = x[0] - 1.0;
        fvec[1] = x[1] - 1.0;
        return 0;
    }
};


int main(int argc, char **argv) {
    /// initialize the objective function
    /// 1st argument is the number of functions, 2nd the number of variables
    MyObjective obj(2, 2);

    /// create an instance of the Levenberg-Marquardt (LM for short) optimizer
    Optimizer_LM lm;

    /// initialized the variables. Later x will be modified after optimization.
    std::vector<double> x = {4.0, -4.0}; // Let's guess the initial values to be (4.0, 4.0)

    /// optimize (i.e., minimizing the objective function).
    bool status = lm.optimize(&obj, x);

    /// retrieve the result.
    std::cout << "the solution is:     " << x[0] << ", " << x[1] << std::endl;
    std::cout << "the expected result: 1, 1" << std::endl;

    return status;
}
