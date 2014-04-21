// ============================================================================= //
// ========================= SETTINGS FILE FOR main.cpp ======================== //
// ============================================================================= //

// ---------------------------------------------------------------------------- //
namespace primal
{
    class InflowBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return x[0] < (0 + DOLFIN_EPS) && on_boundary;
        }
    };

    class GammaSD : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return (x[1] < (0 + DOLFIN_EPS) || x[1] > (2 - DOLFIN_EPS)) && on_boundary;
        }
    };
    
    class NoSlipBoundary : public dolfin::SubDomain
    {
        public:
            NoSlipBoundary () : 
                center (2),
                radius (0.5)
            {
                center[0] = 2.5;
                center[1] = 1;
            }
            
            bool inside (const dolfin::Array<double>& x, bool on_boundary) const
            {
                double dx = x[0] - center[0];
                double dy = x[1] - center[1];
                double r = sqrt (dx * dx + dy * dy);
                
                return r < (radius + 1e-3) && on_boundary;
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
            DirichletBoundary () : 
                center (2),
                radius (0.5)
            {
                center[0] = 2.5;
                center[1] = 1;
            }

            bool inside (const dolfin::Array<double>& x, bool on_boundary) const
            {
                bool onLeftSide  = (x[0] < (0 + DOLFIN_EPS)) && on_boundary;   
                bool onUpperSide = (x[1] > (2 - DOLFIN_EPS)) && on_boundary;   
                bool onLowerSide = (x[1] < (0 + DOLFIN_EPS)) && on_boundary;   

                double dx = x[0] - center[0];
                double dy = x[1] - center[1];
                double r = sqrt (dx * dx + dy * dy);
                bool onCircleBoundary = r < (radius + 1e-3) && on_boundary;

                return onLeftSide || onUpperSide || onLowerSide || onCircleBoundary;
            }
        
        private:
            dolfin::Array<double> center;
            double radius;
    };

    class RobinBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return x[0] > (5 - DOLFIN_EPS) && on_boundary;
        }
    };
    
    class ExternalLoadDomain : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
//            return x[0] >= 1.75 && x[0] <= 3.25 && x[1] >= 0 && x[1] <= 1.75; 
            return x[0] >= 1.75 && x[0] <= 5 && x[1] >= 0 && x[1] <= 2.5; 
        }
    };
}

// ---------------------------------------------------------------------------- //
namespace objective_functional
{
    double sigma_1 = 0.1;
    double sigma_2 = 1.0;
    
    class ControlDomain : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return x[0] < (0 + DOLFIN_EPS) && on_boundary;
        }
    };
    
    class Gradient : public controlproblem::VariableExpression
    {
        public:
        void eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const
        {
            dolfin::Array<double> theta (1);
            evaluateVariable ("theta", theta, x);

            dolfin::Array<double> g (2);
            evaluateVariable ("g", g, x);

            values[0] = sigma_1 * g[0] - theta[0];
            values[1] = sigma_1 * g[1];
        }
        
        std::size_t value_rank () const
        {
            return 1;
        }
        
        std::size_t value_dimension (std::size_t i) const
        {
            return 2;
        }
    };
}

// ---------------------------------------------------------------------------- //

