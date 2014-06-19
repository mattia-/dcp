/* 
 *  Copyright (C) 2014, Mattia Tamellini, mattia.tamellini@gmail.com
 * 
 *  This file is part of the DCP library
 *   
 *   The DCP library is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   The DCP library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with the DCP library.  If not, see <http://www.gnu.org/licenses/>. 
 */ 

// ============================================================================= //
// ========================= SETTINGS FILE FOR main.cpp ======================== //
// ============================================================================= //
// 
#include "function_derivative.h"
#include "diffusion_reaction.h"

class ControlDirichletBC : public dolfin::Expression
{
    public:
        ControlDirichletBC (const dolfin::Function& g) : 
            g_ (g)
        {  }
        
        void eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const
        {
            dolfin::Array <double> x_ (1);
            x_ [0] = x[1];
            g_.eval (values, x_);
        }
        
    private:
        const dolfin::Function& g_;
};

class ValueUpdater
{
    public:
        void operator() (DCP::CompositeDifferentialProblem& compositeProblem, 
                         const dolfin::GenericFunction& dirichletBCValue) const
        {

        }
};


// ---------------------------------------------------------------------------- //
namespace primal
{
    class InflowBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return x[0] <= (0 + DOLFIN_EPS) && on_boundary;
        }
    };

    class GammaSD : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return (x[1] <= (0 + DOLFIN_EPS) || x[1] >= (7 - DOLFIN_EPS)) && on_boundary;
        }
    };
    
    class NoSlipBoundary : public dolfin::SubDomain
    {
        public:
            NoSlipBoundary () : 
                center (2),
                radius (0.5)
            {
                center[0] = 3.5;
                center[1] = 3.5;
            }
            
            bool inside (const dolfin::Array<double>& x, bool on_boundary) const
            {
                double dx = x[0] - center[0];
                double dy = x[1] - center[1];
                double r = sqrt (dx * dx + dy * dy);
                
                return r <= (radius + 1e-3) && on_boundary;
            }

        private:
            dolfin::Array<double> center;
            double radius;
    }; 
}

// ---------------------------------------------------------------------------- //
namespace adjoint
{
    class DirichletBoundary : public dolfin::SubDomain
    {
        public:
            bool inside (const dolfin::Array<double>& x, bool on_boundary) const
            {
                bool onLeftSide  = (x[0] <= (0 + DOLFIN_EPS)) && on_boundary;   
                bool onUpperSide = (x[1] >= (7 - DOLFIN_EPS)) && on_boundary;   
                bool onLowerSide = (x[1] <= (0 + DOLFIN_EPS)) && on_boundary;   

                bool onCircle = circle.inside (x, on_boundary);

                return onLeftSide || onUpperSide || onLowerSide || onCircle;
            }
        
        private:
            primal::NoSlipBoundary circle;
    };

    class RobinBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return x[0] >= (10 - DOLFIN_EPS) && on_boundary;
        }
    };
    
    class ExternalLoadDomain : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
//            return x[0] >= 2.5 && x[0] <= 4.5 && x[1] >= 2.5 && x[1] <= 4.5; 
            return x[0] >= 2.5 && x[0] <= 5 && x[1] >= 2 && x[1] <= 5; 
        }
    };
}

// ---------------------------------------------------------------------------- //
namespace objective_functional
{
    double sigma_1 = 0.01;
    double sigma_2 = 0.01;
    
    class ControlDomain : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return x[0] <= (0 + DOLFIN_EPS) && on_boundary;
        }
    };
    
    class Gradient : public DCP::VariableExpression
    {
        void eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const
        {
            dolfin::Array<double> xFlipped (2);
            xFlipped [0] = x [1];
            xFlipped [1] = x [0];
            dolfin::Array<double> thetaValues (1);
            evaluateVariable ("theta", thetaValues, xFlipped);

            dolfin::Array<double> gValues (1);
            evaluateVariable ("g", gValues, x);

            // compute laplacian of the control g
            dolfin::IntervalMesh mesh (100, 0, 7);
            function_derivative::FunctionSpace V (mesh);
            std::shared_ptr <dolfin::GenericFunction> g = 
                std::const_pointer_cast <dolfin::GenericFunction> (variables_.find ("g") -> second);
            function_derivative::BilinearForm a (V, V);
            function_derivative::LinearForm L (V);
            dolfin::Function gradient (V);
            dolfin::Function laplacian (V);

            int logLevelBackup = dolfin::get_log_level ();
            dolfin::set_log_level (dolfin::WARNING);
            L.set_coefficient ("f", g);

            dolfin::solve (a == L, gradient);

            L.set_coefficient ("f", dolfin::reference_to_no_delete_pointer (gradient));

            dolfin::solve (a == L, laplacian);
            dolfin::set_log_level (logLevelBackup);

            dolfin::Array<double> laplacianValues (1);
            laplacian.eval (laplacianValues, x);

            values[0] = sigma_1 * gValues [0] - sigma_2 * laplacianValues [0] - thetaValues [0];
        }
    };
    
    class SearchDirectionComputer
    {
        public:
            SearchDirectionComputer (const dolfin::Mesh& mesh) :
                mesh_ (mesh),
                functionSpace_ (mesh_),
                a_ (functionSpace_, functionSpace_),
                L_ (functionSpace_)
            { }
            
            void operator() (dolfin::Function& searchDirection, const dolfin::Function& gradient)
            {
                dolfin::Function solution (functionSpace_);
                dolfin::log (dolfin::DBG, "Computing search direction via diffusion-reaction problem...");
                L_.set_coefficient ("f", dolfin::reference_to_no_delete_pointer (gradient));
                dolfin::solve (a_ == L_, solution);
                searchDirection.interpolate (solution);
            }
            
        private:
            const dolfin::Mesh& mesh_;
            diffusion_reaction::FunctionSpace functionSpace_;
            diffusion_reaction::BilinearForm a_;
            diffusion_reaction::LinearForm L_;
    };
}

// ---------------------------------------------------------------------------- //
