#include <dolfin.h>
#include "prova.h"

class RightBd : public dolfin::SubDomain
{
  public :
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return on_boundary
               &&
               ( dolfin::near (x[0],1) );
    }
};
int main (int argc, char * argv[])
{
    dolfin::RectangleMesh mesh (0,0, 1,1, 10,10);
		dolfin::FacetFunction<std::size_t> meshFacets (mesh);
	  meshFacets.set_all(0);
    RightBd noSlipBoundary;
	  noSlipBoundary.mark (meshFacets, 3);
    prova::Form_sigmaLength forma (mesh);
    prova::Form_sigmaLengthR formaR (mesh);
    dolfin::Constant coeff (5.0);
    
    forma.set_coefficient ("coeff", dolfin::reference_to_no_delete_pointer (coeff)); 
    formaR.set_coefficient ("coeff", dolfin::reference_to_no_delete_pointer (coeff)); 
    formaR.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer (meshFacets));
    formaR.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer (meshFacets));
    double val = dolfin::assemble(forma);
    double valR = dolfin::assemble(formaR);
    std::cout << "valore " << val << std::endl;
    std::cout << "valoreR " << valR << std::endl;

return 0;
}
