#ifndef HH__SUBDOMAINTYPE__HH
#define HH__SUBDOMAINTYPE__HH

namespace control_problem
{
	//! \class SubdomainType 
	/*!  \brief Enumeration class used to identify domain parts
	 *
	 *  This class will be used to check which measure to use when setting
	 *  integration subdomains on the linear form and the bilinear form.
	 *  Possible values are:
	 *  TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO   
	 *  \item INTERNAL_CELLS: sets integration subdomains of dimension equal to that of the mesh
	 *  \item BOUNDARY_FACETS: sets integration subdomains of dimension mesh_dimension - 1 on the
	 *  boundary of the mesh
	 *  \item INTERNAL_FACETS: sets integration subdomains of dimension mesh_dimension - 1 inside
	 *  the mesh
	 */
	enum class SubdomainType {INTERNAL_CELLS, BOUNDARY_FACETS, INTERNAL_FACETS};
}
				
#endif
