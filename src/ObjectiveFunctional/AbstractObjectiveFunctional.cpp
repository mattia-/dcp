#include <ObjectiveFunctional/AbstractObjectiveFunctional.hpp>

namespace controlproblem
{
    /************************* CONSTRUCTORS ********************/
    AbstractObjectiveFunctional::AbstractObjectiveFunctional (const boost::shared_ptr <const dolfin::Mesh> mesh) : 
        mesh_ (mesh)
    {
        dolfin::log (dolfin::DBG, "Abstract objective functional created");
    }

    AbstractObjectiveFunctional::AbstractObjectiveFunctional (const dolfin::Mesh& mesh) : 
        mesh_ (new dolfin::Mesh (mesh))
    {
        dolfin::log (dolfin::DBG, "Abstract objective functional created");
    }

    /************************* DESTRUCTOR ********************/
    // this is done for compatibility with gcc 4.6, which doesn't allow virtual members to be defaulted in class body
    AbstractObjectiveFunctional::~AbstractObjectiveFunctional () = default;
}
