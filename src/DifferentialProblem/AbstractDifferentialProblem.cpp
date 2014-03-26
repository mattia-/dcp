#include <DifferentialProblem/AbstractDifferentialProblem.hpp>

namespace control_problem
{
	/************************* CONSTRUCTORS ********************/
	AbstractDifferentialProblem::AbstractDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh,
														  	  const std::shared_ptr<dolfin::FunctionSpace> functionSpace) : 
		mesh_ (mesh),
		functionSpace_ (functionSpace),
		dirichletBCs_ (),
		parameters_ ("differential_problem_parameters"),
		solution_ (*functionSpace)
	{ }



	AbstractDifferentialProblem::AbstractDifferentialProblem (const dolfin::Mesh& mesh,
														  	  const dolfin::FunctionSpace& functionSpace) : 
		mesh_ (new dolfin::Mesh (mesh)),
		functionSpace_ (new dolfin::FunctionSpace (functionSpace)),
		dirichletBCs_ (),
		parameters_ ("differential_problem_parameters"),
		solution_ (functionSpace)
	{ }



	AbstractDifferentialProblem::AbstractDifferentialProblem (dolfin::Mesh&& mesh, 
														  	  dolfin::FunctionSpace&& functionSpace) : 
		mesh_ (std::make_shared<dolfin::Mesh> (mesh)),
		functionSpace_ (std::make_shared<dolfin::FunctionSpace> (functionSpace)),
		dirichletBCs_ (),
		parameters_ ("differential_problem_parameters"),
		solution_ (functionSpace)
	{ }



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



	const dolfin::Parameters& AbstractDifferentialProblem::parameters () const
	{
		return parameters_;
	}



	dolfin::Parameters& AbstractDifferentialProblem::parameters ()
	{
		return parameters_;
	}



	const dolfin::Function& AbstractDifferentialProblem::solution () const
	{
		return solution_;
	}



	/********************** SETTERS ***********************/
	void AbstractDifferentialProblem::addDirichletBC (const dolfin::DirichletBC& dirichletCondition)
	{
		dirichletBCs_.emplace_back (dirichletCondition);
	}



	void AbstractDifferentialProblem::addDirichletBC (dolfin::DirichletBC&& dirichletCondition)
	{
		dirichletBCs_.emplace_back (dirichletCondition);
	}



	void AbstractDifferentialProblem::removeDirichletBC (const std::vector<dolfin::DirichletBC>::iterator& i)
	{
		dirichletBCs_.erase (i);
	}
}
