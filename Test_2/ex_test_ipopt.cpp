#include <iostream>
#include <vector>

#include <ifopt/ipopt_solver.h>
#include <ifopt/problem.h>

#include "test_vars_constr_cost.h"

using namespace ifopt;

#define Num_cost 2

int main(){

    // return type is CppAD::AD<double> actually
    auto constraint_1 = [](const ad_vector &x){
        return -1*x[0]+1*x[1];  
    };

    auto constraint_2 = [](const ad_vector &x){
        return 1*x[0]+(1)*x[1];  
    };

    auto constraint_3 = [](const ad_vector &x){
        return (-1)*x[0]+(-1)*x[1];  
    };

    auto constraint_4 = [](const ad_vector &x){
        return 1*x[0]-3*x[1];  
    };

    auto constraint_5 = [](const ad_vector &x){
        return (x[2]-3)*(x[2]-3) + x[4] -4;
    };

    auto constraint_6 = [](const ad_vector &x){
        return -(x[4]-3)*(x[4]-3) - x[5] + 4;
    };

    auto cost_1 = [](const ad_vector &x){
        return -(25*(x[0]-2)*(x[0]-2) + (x[1]-2)*(x[1]-2) + (x[2]-1)*(x[2]-1) + (x[3]-4)*(x[3]-4) + (x[4]-1)*(x[4]-1));
    };

    auto cost_2 = [](const ad_vector &x){
        return x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[4]*x[4] + x[5]*x[5];
    };
    

#if 0
    //ad_vector x = {1.67, 1.17, 0.17};
    ad_vector x = {1.16, 1.23, 0.60};

    std::cout << "Cost 1: " << cost_1(x) << std::endl;
    std::cout << "Cost 2: " << cost_2(x) << std::endl;
    std::cout << "Cost 3: " << cost_3(x) << std::endl;
#endif

#if 1
    // 定义函数指针类型
    typedef CppAD::AD<double> (*cost_function)(const ad_vector &);

    // 创建函数指针数组
    cost_function cost_functions[] = {cost_1, cost_2};

    Eigen::VectorXd x(Num_var);
    std::vector<Eigen::VectorXd> cost_opt_list;
    std::vector<int> iteration_list;
    
    for(int i = 0; i < Num_cost; i++){
        Problem nlp;

        std::string variable_set_name = "var_set1";
        std::string cost_name = "cost " + std::to_string(i);

        if(i==0){
            x << 1.0, 3.0, 3.0, 3.0, 3.0, 3.0;
            //x << 0.5, 0.5, 0.5, 0.5, 0.5, 0.5;
            nlp.AddVariableSet(std::make_shared<ExVariables>(variable_set_name, x));
        }
        else{
            //std::cout << "x: " << x.transpose() << std::endl;
            nlp.AddVariableSet(std::make_shared<ExVariables>(variable_set_name, x));
        }
        nlp.AddConstraintSet(std::make_shared<ExConstraint<decltype(constraint_1)>>("constraint 1", -2.0, 0, constraint_1));
        nlp.AddConstraintSet(std::make_shared<ExConstraint<decltype(constraint_2)>>("constraint 2", 6.0, 0, constraint_2));
        nlp.AddConstraintSet(std::make_shared<ExConstraint<decltype(constraint_3)>>("constraint 3", 2.0, 0, constraint_3));
        nlp.AddConstraintSet(std::make_shared<ExConstraint<decltype(constraint_4)>>("constraint 4", 2.0, 0, constraint_4));
        nlp.AddConstraintSet(std::make_shared<ExConstraint<decltype(constraint_5)>>("constraint 5", 0.0, 0, constraint_5));
        nlp.AddConstraintSet(std::make_shared<ExConstraint<decltype(constraint_6)>>("constraint 6", 0.0, 0, constraint_6));

        if(i >= 1){
            //std::cout << "cost_opt: " << cost_opt << std::endl;
            int len = cost_opt_list.size();

            for(int j = 0; j < len; j++){
                Eigen::VectorXd cost = cost_opt_list[j];
                std::string cost_constraint_name = "cost_constraint " + std::to_string(j);
                nlp.AddConstraintSet(std::make_shared<ExConstraint<decltype(cost_functions[j])>>(cost_constraint_name, cost(0), 2, cost_functions[j]));
            }
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
        
        int iteration = nlp.GetIterationCount()-2;//! The iteration shown in IPOPT (the acutal iteration) is 2 less than the iteration in ifopt
        std::cout << "Iteration: " << iteration << std::endl;

        Eigen::VectorXd cost_opt = nlp.GetCosts().GetValues();
        std::cout << "Cost: " << cost_opt << std::endl;

        std::cout << "---------------------------------------------------" << std::endl;

        cost_opt_list.push_back(cost_opt);
        iteration_list.push_back(iteration);
    }

    for(int i = 0; i < Num_cost; i++){
        std::cout << "Cost " << i << ": " << cost_opt_list[i].transpose() << std::endl;
        std::cout << "Iteration " << i << ": " << iteration_list[i] << std::endl;
    }

#endif
    
    return 0;
}