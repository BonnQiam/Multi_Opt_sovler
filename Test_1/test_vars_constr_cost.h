#include <vector>

#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>
#include <ifopt/variable_set.h>

#include <cppad/cppad.hpp>

#include <functional>

#define Num_var 3

#define Constraint_type int

typedef CPPAD_TESTVECTOR( CppAD::AD<double> ) ad_vector;
typedef CPPAD_TESTVECTOR( double )            d_vector;

namespace ifopt {
    using Eigen::Vector2d; 
    class ExVariables : public VariableSet {
        private:
            std::vector <double> x_;
        public:
            ExVariables(const std::string& name, VectorXd& x) : VariableSet(3, name){
                //x0_ = x(0); x1_ = x(1); x2_ = x(2);
                for(int i = 0; i < Num_var; i++){
                    x_.push_back(x(i));
                }
            }
            
            void SetVariables(const VectorXd& x) override{
                //x0_ = x(0); x1_ = x(1); x2_ = x(2);
                for(int i = 0; i < Num_var; i++){
                    x_[i] = x(i);
                }
            };
            
            VectorXd GetValues() const override { 
                //VectorXd x(3);
                //x << x0_, x1_, x2_;
                
                VectorXd x(Num_var);
                for(int i = 0; i < Num_var; i++){
                    x(i) = x_[i];
                }

                return x;
            };
            
            VecBound GetBounds() const override{
                VecBound bounds(GetRows());
                //bounds.at(0) = NoBound;
                //bounds.at(1) = NoBound;

                for(int i = 0; i < Num_var; i++){
                    //bounds.at(i) = Bounds(-inf, inf);
                    bounds.at(i) = NoBound;
                }

                bounds.at(2) = Bounds(0.0, inf);
                return bounds;
            }
    };
    
    template <typename Func>
    class ExConstraint : public ConstraintSet {
        private:
            double c;
            Constraint_type t;
            mutable CppAD::ADFun<double> fun_;
        public:
            ExConstraint(const std::string& name, 
            double C, Constraint_type T,
            Func constraint
            ) : ConstraintSet(1, name) {
                c = C; t=T;

                // Independent variable vector
                ad_vector x(Num_var);
                CppAD::Independent(x);

                // Dependent variable vector
                ad_vector y(1);
                y[0] = constraint(x);

                // Create ADFun object
                fun_.Dependent(x, y);
            };

            // 使用 ADFun 对象计算 Constraint
            d_vector compute(const d_vector &x_val) const {
                return fun_.Forward(0, x_val);
            }
            // 使用 ADFun 对象计算 jacobian
            d_vector jacobian(const d_vector &x_val) const {
                return fun_.Jacobian(x_val);
            }
            
            //VectorXd GetValues() const override{
            VectorXd GetValues() const override{
                VectorXd g(GetRows());
                VectorXd x = GetVariables()->GetComponent("var_set1")->GetValues();
                //Vector2d x = GetVariables()->GetComponent("var_set1")->GetValues();

                // Covert VectorXd to d_vector
                d_vector x_val(Num_var);
                for(int i = 0; i < Num_var; i++){
                    x_val[i] = x(i);
                }

                g(0) = compute(x_val)[0];

                return g;
            };

            VecBound GetBounds() const override{
                VecBound b(GetRows());
                //b.at(0) = Bounds(1.0, 1.0);
                if(t == 0)
                // '<' type inequality constraint
                b.at(0) = Bounds(-inf, c);
                else if(t == 1)
                // '>' type inequality constraint
                b.at(0) = Bounds(c, inf);
                else
                // '=' type equality constraint
                b.at(0) = Bounds(c, c);

                return b;
            };

            void FillJacobianBlock(std::string var_set, Jacobian& jac_block) const override{
                if (var_set == "var_set1") {
                    VectorXd x = GetVariables()->GetComponent("var_set1")->GetValues();
                    //Vector2d x = GetVariables()->GetComponent("var_set1")->GetValues();


                    // Covert VectorXd to d_vector
                    d_vector x_val(Num_var);
                    for(int i = 0; i < Num_var; i++){
                        x_val[i] = x(i);
                    }

                    d_vector jac = jacobian(x_val);

                    //jac_block.coeffRef(0, 0) = jac[0];
                    //jac_block.coeffRef(0, 1) = jac[1];
                    //jac_block.coeffRef(0, 2) = jac[2];
                    for(int i = 0; i < Num_var; i++){
                        jac_block.coeffRef(0, i) = jac[i];
                    }
                }
            };
    };

    template <typename Func>
    class ExCost : public CostTerm {
        private:
            mutable CppAD::ADFun<double> fun_;
        public:
            ExCost(const std::string& name,
            Func cost
            ) : CostTerm(name) {
                // Independent variable vector
                ad_vector x(Num_var);
                CppAD::Independent(x);

                // Dependent variable vector
                ad_vector y(1);
                y[0] = cost(x);

                // Create ADFun object
                fun_.Dependent(x, y);
            }

            d_vector compute(const d_vector &x_val) const {
                return fun_.Forward(0, x_val);
            }

            d_vector jacobian(const d_vector &x_val) const {
                return fun_.Jacobian(x_val);
            }

            double GetCost() const override{
                //Vector2d x = GetVariables()->GetComponent("var_set1")->GetValues();
                //return compute({x(0), x(1)})[0];
                VectorXd x = GetVariables()->GetComponent("var_set1")->GetValues();

                // Covert VectorXd to d_vector
                d_vector x_val(Num_var);
                for(int i = 0; i < Num_var; i++){
                    x_val[i] = x(i);
                }

                return compute(x_val)[0];
            };

            void FillJacobianBlock(std::string var_set, Jacobian& jac_block) const override{
                if (var_set == "var_set1") {
                    //Vector2d x = GetVariables()->GetComponent("var_set1")->GetValues();
                    //d_vector jac = jacobian({x(0), x(1)});

                    VectorXd x = GetVariables()->GetComponent("var_set1")->GetValues();
                    // Covert VectorXd to d_vector
                    d_vector x_val(Num_var);
                    for(int i = 0; i < Num_var; i++){
                        x_val[i] = x(i);
                    }

                    d_vector jac = jacobian(x_val);

                    //jac_block.coeffRef(0, 0) = jac[0];
                    //jac_block.coeffRef(0, 1) = jac[1];
                    //jac_block.coeffRef(0, 2) = jac[2];
                    for(int i = 0; i < Num_var; i++){
                        jac_block.coeffRef(0, i) = jac[i];
                    }
                }
            }
    };
}