#include <dcp/utils/DotProduct.h>
#include <dcp/utils/dotproductforms.h>
#include <dolfin/math/basic.h>
#include <dolfin/log/dolfin_log.h>

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
                                const dolfin::Mesh& mesh) const
    {
        std::shared_ptr<dolfin::Form> dotProductComputer = getDotProductComputer_ (first, second, mesh);

        dotProductComputer -> set_coefficient (0, dolfin::reference_to_no_delete_pointer (first));
        dotProductComputer -> set_coefficient (1, dolfin::reference_to_no_delete_pointer (second));

        return dolfin::assemble (*dotProductComputer);
    }



    double DotProduct::compute (const dolfin::Function& first,
                                const dolfin::Function& second) const
    {
        return compute (first, second, *(first.function_space () -> mesh ()));
    }



    double DotProduct::compute (const dcp::TimeDependentFunction& left,
                                const dcp::TimeDependentFunction& right) const
    {
        if (left.size () != right.size ())
        {
            dolfin::dolfin_error ("dcp: DotProduct.cpp",
                                  "perform dot product between objects of type dcp::TimeDependentFunction",
                                  "Size mismatch");
        }

        // loop through time steps and compute dot products
        dolfin::begin (dolfin::DBG, "Computing dot product of time dependent functions...");
        double result = 0;
        double previousDotProduct = 0;
        double currentDotProduct = 0;
        double dt = 0;
        bool isFirstIteration = true;
        for (std::size_t i = 0; i < left.size (); ++i)
        {
            if (dolfin::near (left[i].first, right[i].first) == false)
            {
                dolfin::dolfin_error ("dcp: DotProduct.cpp",
                                      "perform dot product between objects of type dcp::TimeDependentFunction",
                                      "Mismatch in time value at position %d; the two values are %f and %f",
                                      i,
                                      left[i].first,
                                      right[i].first);
            }

            if (isFirstIteration)
            {
                isFirstIteration = false;
                previousDotProduct = compute (left[i].second, right[i].second);
            }
            else
            {
                dt = fabs (left[i].first - left[i-1].first);
                currentDotProduct = compute (left[i].second, right[i].second);
                result += 0.5 * (previousDotProduct + currentDotProduct) * dt;

                previousDotProduct = currentDotProduct;
            }
        }
        dolfin::end (); // Computing dot product of time dependent functions...

        return result;

    }



    double DotProduct::norm (const dolfin::GenericFunction& function,
                             const dolfin::Mesh& mesh) const
    {
        return sqrt (compute (function, function, mesh));
    }



    double DotProduct::norm (const dolfin::Function& function) const
    {
        return norm (function, *(function.function_space () -> mesh ()));
    }



    double DotProduct::norm (const dcp::TimeDependentFunction& function) const
    {
        return sqrt (compute (function, function));
    }



    // ---------------------------------------------------------------------------------------------//



    /********************** PROTECTED METHODS ***********************/
    std::shared_ptr<dolfin::Form> DotProduct::getDotProductComputer_ (const dolfin::GenericFunction& first,
                                                                      const dolfin::GenericFunction& second,
                                                                      const dolfin::Mesh& mesh) const
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
                                      "getDotProductComputer_",
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
                    tmp.reset
                        (new dotproductforms::Form_scalar1D_dotProduct (dolfin::reference_to_no_delete_pointer (mesh)));
                    return tmp;
                }
                else
                {
                    dolfin::dolfin_error ("dcp: DotProduct.cpp",
                                          "getDotProductComputer_",
                                          "No form to compute dot products and norms for mesh cell type \"%s\" and function rank %d",
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
                    tmp.reset
                        (new dotproductforms::Form_scalar2D_dotProduct (dolfin::reference_to_no_delete_pointer (mesh)));
                    return tmp;
                }
                else if (rank == 1)
                {
                    dolfin::log (dolfin::DBG, "Selected vector 2D form to compute dot products and norms");
                    std::shared_ptr <dolfin::Form> tmp = nullptr;
                    tmp.reset
                        (new dotproductforms::Form_vector2D_dotProduct (dolfin::reference_to_no_delete_pointer (mesh)));
                    return tmp;
                }
                else
                {
                    dolfin::dolfin_error ("dcp: BacktrackingOptimizer.cpp",
                                          "getDotProductComputer_",
                                          "No form to compute dot products and norms for mesh cell type \"%s\" and function rank %d",
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
                    tmp.reset
                        (new dotproductforms::Form_scalar3D_dotProduct (dolfin::reference_to_no_delete_pointer (mesh)));
                    return tmp;
                }
                else if (rank == 1)
                {
                    dolfin::log (dolfin::DBG, "Selected vector 3D form to compute dot products and norms");
                    std::shared_ptr <dolfin::Form> tmp = nullptr;
                    tmp.reset
                        (new dotproductforms::Form_vector3D_dotProduct (dolfin::reference_to_no_delete_pointer (mesh)));
                    return tmp;
                }
                else
                {
                    dolfin::dolfin_error ("dcp: BacktrackingOptimizer.cpp",
                                          "getDotProductComputer_",
                                          "No form to compute dot products and norms for mesh cell type \"%s\" and function rank %d",
                                          meshCellType.c_str (),
                                          rank);
                }
            }
        }
        return nullptr;
    }
}
