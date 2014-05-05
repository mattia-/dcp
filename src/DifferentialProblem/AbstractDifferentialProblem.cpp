#include <DifferentialProblem/AbstractDifferentialProblem.hpp>
#include <dolfin/log/dolfin_log.h>
#include <map>
#include <string>
#include <utility>

namespace DCP
{
    /************************* CONSTRUCTORS ********************/
    AbstractDifferentialProblem::AbstractDifferentialProblem (const boost::shared_ptr<dolfin::Mesh> mesh,
                                                              const boost::shared_ptr<dolfin::FunctionSpace> functionSpace) : 
        parameters ("differential_problem_parameters"),
        mesh_ (mesh),
        functionSpace_ (functionSpace),
        dirichletBCs_ (),
        solution_ (*functionSpace_),
        dirichletBCsCounter_ (0)
    { 
        dolfin::log (dolfin::DBG, "AbstractDifferentialProblem object created");
    }



    AbstractDifferentialProblem::AbstractDifferentialProblem (const dolfin::Mesh& mesh,
                                                              const dolfin::FunctionSpace& functionSpace) : 
        parameters ("differential_problem_parameters"),
        mesh_ (new dolfin::Mesh (mesh)),
        functionSpace_ (new dolfin::FunctionSpace (functionSpace)),
        dirichletBCs_ (),
        solution_ (*functionSpace_),
        dirichletBCsCounter_ (0)
    { 
        dolfin::log (dolfin::DBG, "AbstractDifferentialProblem object created"); 
    }



    AbstractDifferentialProblem::AbstractDifferentialProblem (dolfin::Mesh&& mesh, 
                                                              dolfin::FunctionSpace&& functionSpace) : 
        parameters ("differential_problem_parameters"),
        mesh_ (new dolfin::Mesh (std::move (mesh))),
        functionSpace_ (new dolfin::FunctionSpace (std::move (functionSpace))),
        dirichletBCs_ (),
        solution_ (*functionSpace_),
        dirichletBCsCounter_ (0)
    { 
        dolfin::log (dolfin::DBG, "AbstractDifferentialProblem object created"); 
    }



    /********************** GETTERS ***********************/
    const dolfin::Mesh& AbstractDifferentialProblem::mesh () const
    {
        return *mesh_;      
    }



    const dolfin::FunctionSpace& AbstractDifferentialProblem::functionSpace () const
    {
        return *functionSpace_;
    }



    const dolfin::DirichletBC& AbstractDifferentialProblem::dirichletBC (const std::string& bcName) const
    {
        auto bcIterator = dirichletBCs_.find (bcName);
        if (bcIterator == dirichletBCs_.end ())
        {
            dolfin::error ("Cannot find dirichletBC with name \"%s\" in map", bcName.c_str ());
        }
        return bcIterator -> second;
    }



    const std::map<std::string, dolfin::DirichletBC>& AbstractDifferentialProblem::dirichletBCs () const
    {
        return dirichletBCs_;
    }



    const dolfin::Function& AbstractDifferentialProblem::solution () const
    {
        return solution_;
    }



    /********************** SETTERS ***********************/
    bool AbstractDifferentialProblem::addDirichletBC (const dolfin::DirichletBC& dirichletCondition, 
                                                      std::string bcName)
    {
        if (bcName.empty ())
        {
            bcName = "dirichlet_condition_" + std::to_string (dirichletBCsCounter_);
            dirichletBCsCounter_++;
        }
        
        dolfin::log (dolfin::DBG, "Adding dirichlet boundary condition to boundary conditions map with name \"%s\"...",
                     bcName.c_str ());
        auto result = dirichletBCs_.insert (std::make_pair (bcName, dirichletCondition));
        
        if (result.second == false)
        {
            dolfin::warning ("DirichletBC object not inserted because key \"%s\" already in map", bcName.c_str ());
        }
        
        return result.second;
    }



    bool AbstractDifferentialProblem::addDirichletBC (dolfin::DirichletBC&& dirichletCondition,
                                                      std::string bcName)
    {
        std::string bcName_ (bcName);
        if (bcName_.empty ())
        {
            bcName = "dirichlet_condition_" + std::to_string (dirichletBCsCounter_);
            dirichletBCsCounter_++;
        }
        
        dolfin::log (dolfin::DBG, "Adding dirichlet boundary condition to boundary conditions map with name \"%s\"...",
                     bcName.c_str ());
        auto result = dirichletBCs_.insert (std::make_pair (bcName, dirichletCondition));
        
        if (result.second == false)
        {
            dolfin::warning ("DirichletBC object not inserted because key \"%s\" already in map", bcName.c_str ());
        }
        
        return result.second;
    }



    bool AbstractDifferentialProblem::removeDirichletBC (const std::string& bcName)
    {
        dolfin::log (dolfin::DBG, "Removing dirichlet boundary condition \"%s\" from boundary conditions map...", 
                     bcName.c_str ());
        std::size_t nErasedElements = dirichletBCs_.erase (bcName);
        
        if (nErasedElements == 0)
        {
            dolfin::warning ("Dirichlet boundary condition \"%s\" not found in map", bcName.c_str ());
        }
        
        return nErasedElements == 1? true : false;
    }
    


    void AbstractDifferentialProblem::update ()
    {

    }
}
