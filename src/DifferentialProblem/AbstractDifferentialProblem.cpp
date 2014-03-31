#include <DifferentialProblem/AbstractDifferentialProblem.hpp>

namespace control_problem
{
    /************************* CONSTRUCTORS ********************/
    AbstractDifferentialProblem::AbstractDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh,
                                                              const std::shared_ptr<dolfin::FunctionSpace> functionSpace) : 
        parameters ("differential_problem_parameters"),
        mesh_ (mesh),
        functionSpace_ (functionSpace),
        dirichletBCs_ (),
        solution_ (*functionSpace_)
    { 
        dolfin::log (dolfin::DBG, "Abstract differential problem created");
    }



    AbstractDifferentialProblem::AbstractDifferentialProblem (const dolfin::Mesh& mesh,
                                                              const dolfin::FunctionSpace& functionSpace) : 
        parameters ("differential_problem_parameters"),
        mesh_ (new dolfin::Mesh (mesh)),
        functionSpace_ (new dolfin::FunctionSpace (functionSpace)),
        dirichletBCs_ (),
        solution_ (*functionSpace_)
    { 
        dolfin::log (dolfin::DBG, "Abstract differential problem created");
    }



    AbstractDifferentialProblem::AbstractDifferentialProblem (dolfin::Mesh&& mesh, 
                                                              dolfin::FunctionSpace&& functionSpace) : 
        parameters ("differential_problem_parameters"),
        mesh_ (std::make_shared<dolfin::Mesh> (mesh)),
        functionSpace_ (std::make_shared<dolfin::FunctionSpace> (functionSpace)),
        dirichletBCs_ (),
        solution_ (*functionSpace_)
    { 
        dolfin::log (dolfin::DBG, "Abstract differential problem created");
    }



    /************************* DESTRUCTOR ********************/
    // this is done for compatibility with gcc 4.6, which doesn't allow virtual members to be defualted in class body
    AbstractDifferentialProblem::~AbstractDifferentialProblem () = default;



    /********************** GETTERS ***********************/
    const dolfin::Mesh& AbstractDifferentialProblem::mesh () const
    {
        return *mesh_;      
    }



    const dolfin::FunctionSpace& AbstractDifferentialProblem::functionSpace () const
    {
        return *functionSpace_;
    }



    const dolfin::DirichletBC& AbstractDifferentialProblem::dirichletBC (const std::size_t& i) const
    {
        return dirichletBCs_[i];
    }



    const std::vector<dolfin::DirichletBC>& AbstractDifferentialProblem::dirichletBCs () const
    {
        return dirichletBCs_;
    }



    const dolfin::Function& AbstractDifferentialProblem::solution () const
    {
        return solution_;
    }



    /********************** SETTERS ***********************/
    void AbstractDifferentialProblem::addDirichletBC (const dolfin::DirichletBC& dirichletCondition)
    {
        dolfin::log (dolfin::DBG, "Adding dirichlet boundary condition to boundary conditions vector...");
        dirichletBCs_.emplace_back (dirichletCondition);
    }



    void AbstractDifferentialProblem::addDirichletBC (dolfin::DirichletBC&& dirichletCondition)
    {
        dolfin::log (dolfin::DBG, "Adding dirichlet boundary condition to boundary conditions vector...");
        dirichletBCs_.emplace_back (dirichletCondition);
    }



    void AbstractDifferentialProblem::removeDirichletBC (const std::vector<dolfin::DirichletBC>::iterator& i)
    {
        dolfin::log (dolfin::DBG, "Removing dirichlet boundary condition from boundary conditions vector...");
        dirichletBCs_.erase (i);
    }
    


    void AbstractDifferentialProblem::update ()
    {

    }
}
