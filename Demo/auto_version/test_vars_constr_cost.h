#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>
#include <ifopt/variable_set.h>

#include <cppad/cppad.hpp>

#include <functional>

#define Constraint_type int

typedef CPPAD_TESTVECTOR( CppAD::AD<double> ) ad_vector;
typedef CPPAD_TESTVECTOR( double )            d_vector;

namespace ifopt {
    using Eigen::Vector2d; 
    class ExVariables : public VariableSet {
        private:
            double x0_, x1_;
        public:
            ExVariables(const std::string& name, VectorXd& x) : VariableSet(2, name){
                //x0_ = 3.5; x1_ = 1.5;
                x0_ = x(0); x1_ = x(1); 
            }
            
            void SetVariables(const VectorXd& x) override{
                x0_ = x(0); x1_ = x(1);
            };
            
            VectorXd GetValues() const override { 
                return Vector2d(x0_, x1_); 
            };
            
            VecBound GetBounds() const override{
                VecBound bounds(GetRows());
                bounds.at(0) = NoBound;
                bounds.at(1) = NoBound;
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
                ad_vector x(2);
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
                Vector2d x = GetVariables()->GetComponent("var_set1")->GetValues();

                g(0) = compute({x(0), x(1)})[0];

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
                    Vector2d x = GetVariables()->GetComponent("var_set1")->GetValues();

                    d_vector jac = jacobian({x(0), x(1)});

                    jac_block.coeffRef(0, 0) = jac[0];
                    jac_block.coeffRef(0, 1) = jac[1];
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
                ad_vector x(2);
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
                Vector2d x = GetVariables()->GetComponent("var_set1")->GetValues();
                return compute({x(0), x(1)})[0];
            };

            void FillJacobianBlock(std::string var_set, Jacobian& jac_block) const override{
                if (var_set == "var_set1") {
                    Vector2d x = GetVariables()->GetComponent("var_set1")->GetValues();

                    d_vector jac = jacobian({x(0), x(1)});

                    jac_block.coeffRef(0, 0) = jac[0];
                    jac_block.coeffRef(0, 1) = jac[1];
                }
            }
    };
}