#include <iostream>
#include <map>

#include <ifopt/ipopt_solver.h>
#include <ifopt/problem.h>

#include "test_vars_constr_cost.h"

using namespace ifopt;

int main(){

    // return type is CppAD::AD<double> actually
    auto constraint_1 = [](const ad_vector &x){
        return 2*x[0]+x[1];  
    };

    auto constraint_2 = [](const ad_vector &x){
        return 2*x[0]+3*x[1];  
    };

    auto constraint_3 = [](const ad_vector &x){
        return 4*x[0]+3*x[1];  
    };

    auto constraint_4 = [](const ad_vector &x){
        return x[0]+2*x[1];  
    };

    auto cost_1 = [](const ad_vector &x){
        return -8*x[0]-12*x[1];
    };

    auto cost_2 = [](const ad_vector &x){
        return -14*x[0]-10*x[1];  
    };

    // 定义函数指针类型
    typedef CppAD::AD<double> (*cost_function)(const ad_vector &);

    // 创建函数指针数组
    cost_function cost_functions[] = {cost_1, cost_2};

    Eigen::VectorXd x;
    Eigen::VectorXd cost_opt;
    
    for(int i = 0; i < 2; i++){
        Problem nlp;

        std::string variable_set_name = "var_set1";
        std::string cost_constraint_name = "cost_constraint " + std::to_string(i-1);
        std::string cost_name = "cost " + std::to_string(i);

        if(i==0){
            x = Vector2d(3.5, 1.5);
            nlp.AddVariableSet(std::make_shared<ExVariables>(variable_set_name, x));
        }
        else{
            //std::cout << "x: " << x.transpose() << std::endl;
            nlp.AddVariableSet(std::make_shared<ExVariables>(variable_set_name, x));
        }
        nlp.AddConstraintSet(std::make_shared<ExConstraint<decltype(constraint_1)>>("constraint 1", 120.0, 0, constraint_1));
        nlp.AddConstraintSet(std::make_shared<ExConstraint<decltype(constraint_2)>>("constraint 2", 210.0, 0, constraint_2));
        nlp.AddConstraintSet(std::make_shared<ExConstraint<decltype(constraint_3)>>("constraint 3", 270.0, 0, constraint_3));
        nlp.AddConstraintSet(std::make_shared<ExConstraint<decltype(constraint_4)>>("constraint 4", 60.0, 1, constraint_4));

        if(i >= 1){
            //std::cout << "cost_opt: " << cost_opt << std::endl;
            nlp.AddConstraintSet(std::make_shared<ExConstraint<decltype(cost_functions[i-1])>>(cost_constraint_name, cost_opt(0), 2, cost_functions[i-1]));
        }

        nlp.AddCostSet(std::make_shared<ExCost<decltype(cost_functions[i])>>(cost_name, cost_functions[i]));
        
        nlp.PrintCurrent();

        // 2. choose solver and options
        IpoptSolver ipopt;
        ipopt.SetOption("linear_solver", "ma27");
        ipopt.SetOption("jacobian_approximation", "exact");

        // 3 . solve
        ipopt.Solve(nlp);
        
        x = nlp.GetOptVariables()->GetValues();
        std::cout << x.transpose() << std::endl;

        int iteration = nlp.GetIterationCount();
        std::cout << "Iteration: " << iteration << std::endl;

        cost_opt = nlp.GetCosts().GetValues();
        std::cout << "Cost: " << cost_opt << std::endl;
    }

    return 0;
}