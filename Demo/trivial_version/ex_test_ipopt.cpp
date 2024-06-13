/******************************************************************************
Copyright (c) 2017, Alexander W. Winkler, ETH Zurich. All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the name of ETH ZURICH nor the names of its contributors may be
      used to endorse or promote products derived from this software without
      specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************/

#include <iostream>

#include <ifopt/ipopt_solver.h>
#include <ifopt/problem.h>

//#include <ifopt/test_vars_constr_cost.h>

//#include <ipopt_solver.h>
//#include <problem.h>

#include "test_vars_constr_cost.h"
//#include "test_vars_constr_cost_bak.h"


using namespace ifopt;

int main()
{
  // 1. define the problem
  Problem nlp;
  nlp.AddVariableSet(std::make_shared<ExVariables>());
  // std::make_shared<ExVariables>() means create a shared pointer to a new ExVariables object

  nlp.AddConstraintSet(std::make_shared<ExConstraint>("constraint 1", 2.0, 1.0, 120.0, 0));
  nlp.AddConstraintSet(std::make_shared<ExConstraint>("constraint 2", 2.0, 3.0, 210.0, 0));
  nlp.AddConstraintSet(std::make_shared<ExConstraint>("constraint 3", 4.0, 3.0, 270.0, 0));
  nlp.AddConstraintSet(std::make_shared<ExConstraint>("constraint 4", 1.0, 2.0, 60.0, 1));

  //nlp.AddConstraintSet(std::make_shared<ExConstraint>("constraint 5", 8.0, 12.0, 840.0, 2));
  nlp.AddConstraintSet(std::make_shared<ExConstraint>("constraint 5", -8.0, -12.0, -840.0, 2));

  nlp.AddCostSet(std::make_shared<ExCost>());
  
  nlp.PrintCurrent();

  // 2. choose solver and options
  IpoptSolver ipopt;
  ipopt.SetOption("linear_solver", "ma27");
  ipopt.SetOption("jacobian_approximation", "exact");

  // 3 . solve
  ipopt.Solve(nlp);
  
  Eigen::VectorXd x = nlp.GetOptVariables()->GetValues();
  std::cout << x.transpose() << std::endl;

  int iteration = nlp.GetIterationCount();
  std::cout << "Iteration: " << iteration << std::endl;

  Eigen::VectorXd cost_opt = nlp.GetCosts().GetValues();
  std::cout << "Cost: " << cost_opt << std::endl;

  /*
  // 4. test if solution correct
  double eps = 1e-5;  //double precision
  assert(1.0 - eps < x(0) && x(0) < 1.0 + eps); // This meangs that x(0) should be 1.0
  assert(0.0 - eps < x(1) && x(1) < 0.0 + eps);
  */
}
