#include <dcp/utils/DotProduct.h>
#include <dcp/utils/dotproductforms.h>

namespace dcp
{
    /************************* CONSTRUCTORS ********************/
    DotProduct::DotProduct () : 
        dotProductComputer_ (nullptr)
    {
        dolfin::log (dolfin::DBG, "DotProduct object created");
    };
    

    
    /********************** METHODS ***********************/
    void DotProduct::setDotProductComputer (const dolfin::Form& dotProductComputer)
    {
        dotProductComputer_.reset (new dolfin::Form (dotProductComputer));
    }
    


    void DotProduct::resetDotProductComputer ()
    {
        dotProductComputer_.reset ();
    }
    


    double DotProduct::compute (const dolfin::GenericFunction& first, 
                                const dolfin::GenericFunction& second,
                                const dolfin::Mesh& mesh)
    {
        std::shared_ptr<dolfin::Form> dotProductComputer = getDotProductComputer (first, second, mesh);
        
        dotProductComputer -> set_coefficient (0, dolfin::reference_to_no_delete_pointer (first));
        dotProductComputer -> set_coefficient (1, dolfin::reference_to_no_delete_pointer (second));
        
        return dolfin::assemble (*dotProductComputer);
    }
    


    double DotProduct::norm (const dolfin::GenericFunction& function, 
                             const dolfin::Mesh& mesh)
    {
        return sqrt (compute (function, function, mesh));
    }
    


    // ---------------------------------------------------------------------------------------------//



    std::shared_ptr<dolfin::Form> DotProduct::getDotProductComputer (const dolfin::GenericFunction& first,
                                                                     const dolfin::GenericFunction& second,
                                                                     const dolfin::Mesh& mesh)
    {
        if (dotProductComputer_ != nullptr)
        {
            dolfin::log (dolfin::DBG, "Using externally set protected member variable to compute dot products and norms");
            return dotProductComputer_;
        }
        else
        {
            // get type of cells in mesh
            std::string meshCellType = dolfin::CellType::type2string (mesh.type ().cell_type ());
            dolfin::log (dolfin::DBG, "Mesh cell type is: %s", meshCellType.c_str ());

            // get rank of first function
            int firstFunctionRank = first.value_rank ();
            dolfin::log (dolfin::DBG, "First function rank is: %d", firstFunctionRank);
            
            // get rank of second function
            int secondFunctionRank = second.value_rank ();
            dolfin::log (dolfin::DBG, "Second function rank is: %d", secondFunctionRank);
            
            // check if ranks match
            if (firstFunctionRank != secondFunctionRank)
            {
                dolfin::dolfin_error ("dcp: DotProduct.cpp",
                                      "getDotProductComputer",
                                      "Ranks of first and second function do not match");
                return nullptr;
            }
            
            int rank = firstFunctionRank; // it's the same as secondFunctionRank anyway
            
            if (meshCellType == "interval")
            {
                if (rank == 0)
                {
                    dolfin::log (dolfin::DBG, "Selected scalar 1D form to compute dot products and norms");
                    std::shared_ptr <dolfin::Form> tmp = nullptr;
                    tmp.reset (new dotproductforms::Form_scalar1D_dotProduct (mesh));
                    return tmp;
                }
                else
                {
                    dolfin::dolfin_error ("dcp: DotProduct.cpp", 
                                          "getDotProductComputer",
                                          "No form to compute dot products and norms for mesh cell type \"%s\" and", 
                                          "function rank %d", 
                                          meshCellType.c_str (),
                                          rank);
                }
            }

            if (meshCellType == "triangle")
            {
                if (rank == 0)
                {
                    dolfin::log (dolfin::DBG, "Selected scalar 2D form to compute dot products and norms");
                    std::shared_ptr <dolfin::Form> tmp = nullptr;
                    tmp.reset (new dotproductforms::Form_scalar2D_dotProduct (mesh));
                    return tmp;
                }
                else if (rank == 1)
                {
                    dolfin::log (dolfin::DBG, "Selected vector 2D form to compute dot products and norms");
                    std::shared_ptr <dolfin::Form> tmp = nullptr;
                    tmp.reset (new dotproductforms::Form_vector2D_dotProduct (mesh));
                    return tmp;
                }
                else
                {
                    dolfin::dolfin_error ("dcp: BacktrackingOptimizer.cpp", 
                                          "getDotProductComputer",
                                          "No form to compute dot products and norms for mesh cell type \"%s\" and", 
                                          "function rank %d", 
                                          meshCellType.c_str (),
                                          rank);
                }
            }

            if (meshCellType == "tetrahedron")
            {
                if (rank == 0)
                {
                    dolfin::log (dolfin::DBG, "Selected scalar 3D form to compute dot products and norms");
                    std::shared_ptr <dolfin::Form> tmp = nullptr;
                    tmp.reset (new dotproductforms::Form_scalar3D_dotProduct (mesh));
                    return tmp;
                }
                else if (rank == 1)
                {
                    dolfin::log (dolfin::DBG, "Selected vector 3D form to compute dot products and norms");
                    std::shared_ptr <dolfin::Form> tmp = nullptr;
                    tmp.reset (new dotproductforms::Form_vector3D_dotProduct (mesh));
                    return tmp;
                }
                else
                {
                    dolfin::dolfin_error ("dcp: BacktrackingOptimizer.cpp", 
                                          "getDotProductComputer",
                                          "No form to compute dot products and norms for mesh cell type \"%s\" and", 
                                          "function rank %d", 
                                          meshCellType.c_str (),
                                          rank);
                }
            }
        }
        return nullptr;
    }
}
