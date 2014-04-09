#include <ObjectiveFunctional/AbstractObjectiveFunctional.hpp>

namespace controlproblem
{
    /************************* CONSTRUCTORS ********************/
    AbstractObjectiveFunctional::AbstractObjectiveFunctional (const boost::shared_ptr <const dolfin::Mesh> mesh) : 
        mesh_ (mesh)
    {
        dolfin::log (dolfin::DBG, "AbstractObjectiveFunctional object created");
    }

    AbstractObjectiveFunctional::AbstractObjectiveFunctional (const dolfin::Mesh& mesh) : 
        mesh_ (new dolfin::Mesh (mesh))
    {
        dolfin::log (dolfin::DBG, "AbstractObjectiveFunctional object created");
    }
    
    

    /******************* GETTERS *******************/
    const dolfin::Mesh& AbstractObjectiveFunctional::mesh () const
    {
        return *mesh_;
    }

}
