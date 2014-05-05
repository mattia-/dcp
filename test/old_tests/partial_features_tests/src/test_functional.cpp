#include "Poisson.h"
#include "NavierStokes.h"
#include "Functional.h"
#include <memory>
#include "DifferentialProblem/LinearDifferentialProblem.hpp"
#include "DifferentialProblem/NonlinearDifferentialProblem.hpp"
#include "DifferentialProblem/CompositeDifferentialProblem.hpp"
#include "Utils/SubdomainType.hpp"
#include "ObjectiveFunctional/VariableExpression.hpp"
#include "ObjectiveFunctional/ObjectiveFunctional.hpp"
#include <iostream>
#include <dolfin.h>

namespace Poisson
{
    class ExternalLoad : public dolfin::Expression
    {
        void eval (dolfin::Array<double> & values, const dolfin::Array<double> & x) const
        {
            values [0] = 10 * exp (-((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5)) / 0.02);
        }
    };
    
    class NeumannCondition : public dolfin::Expression
    {
        void eval (dolfin::Array<double> & values, const dolfin::Array<double> & x) const
        {
            values [0] = sin (5 * x [0]);
        }
    };
    
    class NeumannBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double> & x, bool on_boundary) const
        {
            return (x[1] < (0 + DOLFIN_EPS) || x[1] > (1 - DOLFIN_EPS)) && on_boundary;
        }
    };
    
    class DirichletCondition : public dolfin::Expression
    {
        void eval (dolfin::Array<double> & values, const dolfin::Array<double> & x) const
        {
            values [0] = 0;
        }
    };

    class DirichletBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double> & x, bool on_boundary) const
        {
            return (x[0] < 0 + DOLFIN_EPS || x[0] > 1 - DOLFIN_EPS) && on_boundary;
        }
    };
    
    class UnitaryConstant : public dolfin::Expression
    {
        void eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const
        {
            values[0] = 1.0;
        }
    };
}

namespace NavierStokes
{
    class MovingLidBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double> & x, bool on_boundary) const
        {
            return x [1] > (1.0 - DOLFIN_EPS) && on_boundary;
        }
    };

    class NoSlipBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double> & x, bool on_boundary) const
        {
            return x [1] < (1.0 - DOLFIN_EPS) && on_boundary;
        }
    };
}

class ProvaVE : public DCP::VariableExpression
{
    private:
        bool inside (const dolfin::Array<double>& x) const
        {
            return x[0] < 1 && x [1] > 0;
        }
    public:
        void eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const
        {
            dolfin::Array<double> bValues (3);
            dolfin::Function bb (function ("b"));
            bb.eval (bValues, x);

            dolfin::Array<double> cValues (3);
            evaluateVariable ("c", cValues, x);
            
            std::cout << "valori a caso:" << std::endl;
            std::cout << value_rank () << std::endl;
            std::cout << inside (x) << std::endl;
            
            dolfin::cout << "bValues: " << dolfin::endl;
            dolfin::cout << bValues[0] << dolfin::endl;
            dolfin::cout << bValues[1] << dolfin::endl;
            dolfin::cout << bValues[2] << dolfin::endl;
            dolfin::cout << dolfin::endl;

            dolfin::cout << "cValues: " << dolfin::endl;
            dolfin::cout << cValues[0] << dolfin::endl;
            dolfin::cout << cValues[1] << dolfin::endl;
            dolfin::cout << cValues[2] << dolfin::endl;
            dolfin::cout << dolfin::endl;
            
            values [0] = 5 + bValues [0] + cValues [0];
            values [1] = 4 * x [0] + x [1] * x [2];
            values [2] = 4 * x [0] + x [1] * x [2];
        }
};


int main ()
{
    dolfin::set_log_level (dolfin::DBG);
    dolfin::UnitSquareMesh mesh (10, 10);
    Poisson::FunctionSpace V (mesh);
    DCP::VariableExpression* a = new ProvaVE;
    dolfin::Array<double> x (3);
    x[0] = 0.1;
    x[1] = 0.2;
    x[2] = 0.3;

    dolfin::Array<double> values (3);

    dolfin::Function b (V);

    a->setCoefficient ("b", dolfin::reference_to_no_delete_pointer (b));
    a->setCoefficient ("c", dolfin::reference_to_no_delete_pointer (b));
    a->setCoefficient ("b", dolfin::reference_to_no_delete_pointer (b));

    a->eval (values, x);

    dolfin::cout << "values: " << dolfin::endl;
    dolfin::cout << values[0] << dolfin::endl;
    dolfin::cout << values[1] << dolfin::endl;
    dolfin::cout << values[2] << dolfin::endl;
    dolfin::cout << dolfin::endl;

    Functional::Functional aaa (mesh);
    aaa.set_coefficient ("f", dolfin::reference_to_no_delete_pointer (b));
    std::cout << dolfin::assemble (aaa) << std::endl;

    // ------ with functional class ------ //
    ProvaVE a2;
    a2.setCoefficient ("b", dolfin::reference_to_no_delete_pointer (b));
    a2.setCoefficient ("c", dolfin::reference_to_no_delete_pointer (b));
    a2.setCoefficient ("b", dolfin::reference_to_no_delete_pointer (b));
    
    dolfin::cout << dolfin::endl;
    dolfin::cout << "PRIMO OGGETTO" << dolfin::endl;
    dolfin::cout << dolfin::endl;
    
    DCP::ObjectiveFunctional<Functional::Functional, ProvaVE> objF (reference_to_no_delete_pointer (mesh));

    objF.setCoefficient ("functional", dolfin::reference_to_no_delete_pointer (b), "f");
    objF.setCoefficient ("gradient", dolfin::reference_to_no_delete_pointer (b), "b");
    
    boost::shared_ptr<dolfin::Constant> con (new dolfin::Constant (42));
    objF.setCoefficient ("gradient", con, "c");

    dolfin::cout << objF.evaluateFunctional () << dolfin::endl;
    
    objF.evaluateGradient (values, x);

    dolfin::cout << "values: " << dolfin::endl;
    dolfin::cout << values[0] << dolfin::endl;
    dolfin::cout << values[1] << dolfin::endl;
    dolfin::cout << values[2] << dolfin::endl;
    dolfin::cout << dolfin::endl;
    
    dolfin::cout << dolfin::endl;
    dolfin::cout << "SECONDO OGGETTO (COPIA DEL PRIMO)" << dolfin::endl;
    dolfin::cout << dolfin::endl;
    
    DCP::ObjectiveFunctional<Functional::Functional, ProvaVE> objF2 (objF);

    dolfin::cout << objF2.evaluateFunctional () << dolfin::endl;
    
    objF2.evaluateGradient (values, x);

    dolfin::cout << "values: " << dolfin::endl;
    dolfin::cout << values[0] << dolfin::endl;
    dolfin::cout << values[1] << dolfin::endl;
    dolfin::cout << values[2] << dolfin::endl;
    dolfin::cout << dolfin::endl;
    
    
    a->setCoefficient ("c", con);
    dolfin::Function dunz1 (V);
    dunz1 = *a;
    Poisson::UnitaryConstant cccc;
    dolfin::Function dunz2 (V);
    dunz2 = cccc;
    dolfin::Function dunz3 (V);
    dunz3 = dunz1 + dunz2;
    dolfin::plot (dunz3);
//    dolfin::plot (dolfin::grad(dunz3));
//    dolfin::plot (dolfin::nabla_grad(dunz3));
    dolfin::interactive ();
    return 0;
}
