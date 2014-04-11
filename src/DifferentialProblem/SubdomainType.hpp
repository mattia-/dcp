#ifndef SRC_DIFFERENTIALPROBLEM_SUBDOMAINTYPE_HPP_INCLUDE_GUARD
#define SRC_DIFFERENTIALPROBLEM_SUBDOMAINTYPE_HPP_INCLUDE_GUARD

namespace controlproblem
{
    /*! \brief Enumeration class used to identify domain parts
     *
     *  This class will be used to check which measure to use when setting
     *  integration subdomains on the linear form and the bilinear form.
     *  Possible values are:
     *  \li INTERNAL_CELLS: sets integration subdomains of dimension equal to that of the mesh
     *  \li BOUNDARY_FACETS: sets integration subdomains of dimension mesh_dimension - 1 on the
     *  boundary of the mesh
     *  \li INTERNAL_FACETS: sets integration subdomains of dimension mesh_dimension - 1 inside
     *  the mesh
     */
    enum class SubdomainType {INTERNAL_CELLS, BOUNDARY_FACETS, INTERNAL_FACETS};
}
                
#endif
