// Begin demo

#include <dolfin.h>
#include "simpleNavierStokes.h"

/*// Initial conditions
class InitialConditions : public dolfin::Expression
{
public:

  InitialConditions() : dolfin::Expression(2)
  {
    dolfin::seed(2 + dolfin::MPI::rank(MPI_COMM_WORLD));
  }

  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
  {
    values[0]= 1.0;//0.63 + 0.02*(0.5 - dolfin::rand());
    values[1]= 0.0;
  }

};*/

// Subdomains and expressions for BCs
class MovingLid : public dolfin::SubDomain
{
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return dolfin::near (x[1], 1) && !(dolfin::near (x[0],0)) && on_boundary;
    }
};

class FixedWalls : public dolfin::SubDomain
{
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return (dolfin::near (x[1], 0) && on_boundary)
               ||
               (dolfin::near (x[0], 0) && on_boundary) 
               ||
               (dolfin::near (x[0], 1) && on_boundary);
               
    }
};
class DirichletBoundary : public dolfin::SubDomain
{
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return  on_boundary
                && (
                dolfin::near(x[1], 0)
                ||
                dolfin::near(x[0], 0)
                ||
                dolfin::near(x[0], 1) );
    }
};
class DirichletData : public dolfin::Expression
{
public:

  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
  {
    values[0]= 0.0;// * x[0] * (1-x[0]);
    values[1]= 0.2 * x[0] * (1-x[0]) * (0.5-x[0]);
  }

  std::size_t value_rank() const
  {
    return 1;
  }

  std::size_t value_dimension(std::size_t i) const
  {
    return 2;
  }

};
class DirichletDataComplete : public dolfin::Expression
{
public:

  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
  {
    values[0]= 0.0;// * x[0] * (1-x[0]);
    values[1]= 0.01 * x[0] * (1-x[0]) * x[1];
    values[2]= 0.0;
  }

  std::size_t value_rank() const
  {
    return 1;
  }

  std::size_t value_dimension(std::size_t i) const
  {
    return 3;
  }

};

// User defined nonlinear problem
template <class T_FunctionSpace, class T_ResidualForm, class T_JacobianForm>
class UflToNewton : public dolfin::NonlinearProblem
{
  public:

    // Constructor
    UflToNewton (const dolfin::Mesh& mesh, const dolfin::Constant& nu,
                 const dolfin::Constant& gamma, const dolfin::Constant& dt) :
        functionSpace_ (new T_FunctionSpace(mesh)),
        F_ (new T_ResidualForm (*functionSpace_)),
        J_ (new T_JacobianForm (*functionSpace_, *functionSpace_)),
        w_ (new dolfin::Function((*functionSpace_)[0]->collapse())),
        solution_ (new dolfin::Function(functionSpace_)),
//notrial        dSolution_ (new dolfin::Function(functionSpace_)),
        previousSolution_ (new dolfin::Function(functionSpace_))
    {
/*      // Initialize class
      // Unfortunately C++ does not allow namespaces as template arguments
      init<simpleNavierStokes::FunctionSpace, simpleNavierStokes::JacobianForm,
           simpleNavierStokes::ResidualForm>(mesh, nu, gamma, dt, w);
*/
        F_ -> set_coefficient ("nu", dolfin::reference_to_no_delete_pointer(nu));
//        F_ -> set_coefficient ("gamma", dolfin::reference_to_no_delete_pointer(gamma));
        F_ -> set_coefficient ("dt", dolfin::reference_to_no_delete_pointer(dt));
//        F_ -> set_coefficient ("w", w_);
        F_ -> set_coefficient ("u_old", dolfin::reference_to_no_delete_pointer((*previousSolution_)[0]));
        F_ -> set_coefficient ("trial", solution_);

        J_ -> set_coefficient ("nu", dolfin::reference_to_no_delete_pointer(nu));
        J_ -> set_coefficient ("dt", dolfin::reference_to_no_delete_pointer(dt));
//        J_ -> set_coefficient ("w", w_);
        J_ -> set_coefficient ("trial", solution_);
//notrial        J_ -> set_coefficient ("dtrial", *dSolution_);

        dolfin::Constant zero(0,0,0);
        DirichletDataComplete initialSolution;
        *solution_ = initialSolution;
        *previousSolution_ = zero;
/*        dolfin::Constant w(0,0);
        *w_ = w;*/
    }

    void addDirichletBC (std::string&& str, dolfin::DirichletBC&& dirichletBC)
    {
        dirichletBCs_.insert (std::make_pair(str,dirichletBC));
    }

    // User defined residual vector
    void F(dolfin::GenericVector& b, const dolfin::GenericVector& x)
    {
      // Assemble RHS (Neumann boundary conditions)
      dolfin::Assembler assembler;
      assembler.assemble (b, *F_);
      for (auto it = dirichletBCs_.begin(); it != dirichletBCs_.end(); it++)
          (it->second).apply (b, x);
    }

    // User defined assemble of Jacobian
    void J(dolfin::GenericMatrix& A, const dolfin::GenericVector& x)
    {
      // Assemble system
      dolfin::Assembler assembler;
      assembler.assemble(A, *J_);
      for (auto it = dirichletBCs_.begin(); it != dirichletBCs_.end(); it++)
          (it->second).apply (A);
    }

    // Return solution function
    dolfin::Function& solution()
    { return *solution_; }

    // Return solution function
/*    Function& u0()
    { return *_u0; }
*/
    dolfin::FunctionSpace& functionSpace()
    { return *functionSpace_; }

  private:

//    template<class T_FunctionSpace, class T_JacobianForm, class T_ResidualForm>
    void init (const dolfin::Mesh& mesh, const dolfin::Constant& nu,
               const dolfin::Constant& gamma, const dolfin::Constant& dt,
               const dolfin::Function& w)
    {
std::cerr << "HAI USATO LA INIT!!!" << std::endl; exit(1);
/*      // Create function space and functions
      std::shared_ptr<T_FunctionSpace> V(new T_FunctionSpace(mesh));
      solution_.reset(new Function(V));
      previousSolution_.reset(new Function(V));

      // Create forms and attach functions
      T_JacobianForm* J_ = new T_JacobianForm(V, V);
      T_ResidualForm* F_ = new T_ResidualForm(V);
      _a->u = *_u;
      _a->lmbda = lambda; _a->dt = dt; _a->theta = theta;
      _L->u = *_u; _L->u0 = *_u0;
      _L->lmbda = lambda; _L->dt = dt; _L->theta = theta;

      // Wrap pointers in a smart pointer
      a.reset(_a);
      L.reset(_L);

      // Set solution to intitial condition
      InitialConditions u_initial;
      *_u = u_initial;
*/
    }

    // Function space, forms and functions
    std::shared_ptr<dolfin::FunctionSpace> functionSpace_;
    std::shared_ptr<dolfin::Form> F_;
    std::shared_ptr<dolfin::Form> J_;
    std::shared_ptr<dolfin::Function> w_;
    std::shared_ptr<dolfin::Function> solution_;
    std::shared_ptr<dolfin::Function> previousSolution_;
    std::map <std::string, dolfin::DirichletBC> dirichletBCs_;
};


int main(int argc, char* argv[])
{
  dolfin::init(argc, argv);

  // Mesh
//  dolfin::UnitSquareMesh mesh(20, 20);
  dolfin::RectangleMesh mesh(0,0, 1,1, 20, 20);

  // Time stepping and model parameters
  dolfin::Constant dt(0.1);
  dolfin::Constant nu(1.0e-01);
  dolfin::Constant gamma(0.0);//7.3e-5);

  double t = 0.0;
  double T = 10*dt;

  // Boundary conditions
  dolfin::Constant movingLidVelocity (1.0, 0.0);
//  DirichletData movingLidVelocity;
  dolfin::Constant noSlipDirichletBC (0.0, 0.0);
  MovingLid movingLid;
  FixedWalls fixedWalls;

  // Create user-defined nonlinear problem
  UflToNewton < simpleNavierStokes::FunctionSpace, simpleNavierStokes::ResidualForm, simpleNavierStokes::JacobianForm >
       nonlinearProblem(mesh, nu, gamma, dt);
  nonlinearProblem.addDirichletBC (std::string("moving"), dolfin::DirichletBC (* nonlinearProblem.functionSpace()[0], movingLidVelocity, movingLid));
  nonlinearProblem.addDirichletBC (std::string("fixed"), dolfin::DirichletBC (* nonlinearProblem.functionSpace()[0], noSlipDirichletBC, fixedWalls));

  // Solution functions
  dolfin::Function& solution = nonlinearProblem.solution();
  dolfin::Function& previousSolution = nonlinearProblem.solution();

  // Create nonlinear solver and set parameters
  dolfin::NewtonSolver newton_solver;
  newton_solver.parameters["linear_solver"] = "lu";
  newton_solver.parameters["convergence_criterion"] = "incremental";
  newton_solver.parameters["maximum_iterations"] = 50;
/*  newton_solver.parameters["relative_tolerance"] = 1e-6;
  newton_solver.parameters["absolute_tolerance"] = 1e-15;*/

  // Save initial condition to file
  dolfin::File file("cahn_hilliard.pvd");
  file << std::make_pair(&solution, t);

  // Solve
  while (t < T)
  {
    // Update for next time step
    t += dt;
    *previousSolution.vector() = *solution.vector();

    // Solve
    newton_solver.solve(nonlinearProblem, *solution.vector());

    // Save function to file
    file << std::pair<const dolfin::Function*, double>(&solution, t);
  }

  // Plot solution
  dolfin::plot(solution[0]);
  dolfin::interactive();

  return 0;
}

