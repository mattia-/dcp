/* 
 *  Copyright (C) 2015, Ivan Fumagalli, ivan.fumagalli@polimi.it
 * 
 *  This file is part of the DCP library
 *   
 *   The DCP library is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   The DCP library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with the DCP library.  If not, see <http://www.gnu.org/licenses/>. 
 */ 

#ifndef IVAN_UFLTONEWTON_H_INCLUDE_GUARD
#define IVAN_UFLTONEWTON_H_INCLUDE_GUARD

#include <dolfin.h>
#include <dcp/differential_problems/SubdomainType.h>
#include "utilities.h"

namespace Ivan
{

  // User defined nonlinear problem
  template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_ResidualForm, class T_JacobianForm, class T_AdditionalForm>
  class UflToNewton : public dolfin::NonlinearProblem
  {

    public:
  
      // Constructor
      UflToNewton (const dolfin::Mesh& mesh) :
          functionSpace_ (new T_FunctionSpace(mesh)),
          additionalFunctionSpace_ (new T_AdditionalFunctionSpace(mesh)),
          residualForm_ (T_ResidualForm (*functionSpace_)),
          jacobianForm_ (T_JacobianForm (*functionSpace_, *functionSpace_)),
          additionalForm_ (T_AdditionalForm (*additionalFunctionSpace_)),
          solution_ (new dolfin::Function(functionSpace_)),
          previousSolution_ (new dolfin::Function(functionSpace_))
      {
  /*      // Initialize class
        // Unfortunately C++ does not allow namespaces as template arguments
        init<simpleNavierStokes::FunctionSpace, simpleNavierStokes::JacobianForm,
             simpleNavierStokes::ResidualForm>(mesh, nu, gamma, dt, w);
  */
      }
  
      void addDirichletBC (std::string&& str, dolfin::DirichletBC&& dirichletBC)
      {
          dirichletBCs_.insert (std::make_pair(str,dirichletBC));
      }
  
      // User defined residual vector
      void F(dolfin::GenericVector& b, const dolfin::GenericVector& x)
      {
        // Assemble RHS (Neumann boundary conditions)
        dolfin::Assembler assembler;
        assembler.assemble (b, residualForm_);
        dolfin::Vector vec(MPI_COMM_WORLD,b.size());
        assembler.assemble (vec, additionalForm_);
std::cerr << ".h OK fino a " << __LINE__ << std::endl;
        processAdditionalVector (vec);
std::cerr << ".h OK fino a " << __LINE__ << std::endl;
for (auto i=0; i!=vec.size(); ++i) std::cerr << vec[i] << ", "; std::cerr << std::endl;
for (auto i=0; i!=b.size(); ++i) std::cerr << b[i] << ", "; std::cerr << std::endl;
        b -= vec;
std::cerr << ".h OK fino a " << __LINE__ << std::endl;
        for (auto it = dirichletBCs_.begin(); it != dirichletBCs_.end(); it++)
            (it->second).apply (b, x);
std::cerr << ".h OK fino a " << __LINE__ << std::endl;
dolfin::Function bFun (*additionalFunctionSpace_);
* bFun.vector() = b;
dolfin::plot (bFun[0], "velocity");
dolfin::plot (bFun[0][0],"x velocity");
dolfin::plot (bFun[0][1],"y velocity"); dolfin::interactive();
std::cerr << std::endl; exit(0);
      }
  
      // User defined assemble of Jacobian
      void J(dolfin::GenericMatrix& A, const dolfin::GenericVector& x)
      {
        // Assemble system
        dolfin::Assembler assembler;
        assembler.assemble(A, jacobianForm_);
        for (auto it = dirichletBCs_.begin(); it != dirichletBCs_.end(); it++)
            (it->second).apply (A);
      }
  
      // Return solution function
      dolfin::Function& solution()
      { return *solution_; }
  
      // Return previous solution function
      dolfin::Function& previousSolution()
      { return *previousSolution_; }
  
      dolfin::FunctionSpace& functionSpace()
      { return *functionSpace_; }
  
      virtual void setCoefficient (const std::string& coefficientType, 
                                   const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                   const std::string& coefficientName = "default");// override;

      virtual void setIntegrationSubdomain (const std::string& formType,
                                             std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                             const dcp::SubdomainType& subdomainType);// override;
                
    private:
  
  //    template<class T_FunctionSpace, class T_JacobianForm, class T_ResidualForm>
      void init (const dolfin::Mesh& mesh, const dolfin::Constant& nu,
                 const dolfin::Constant& gamma, const dolfin::Constant& dt,
                 const dolfin::Function& w)
      {
  std::cerr << "HAI USATO LA INIT!!!" << std::endl; exit(1);
  /*      // Create function space and functions
        std::shared_ptr<T_FunctionSpace> V(new T_FunctionSpace(mesh));
        solution_.reset(new Function(V));
        previousSolution_.reset(new Function(V));
  
        // Create forms and attach functions
        T_JacobianForm* jacobianForm_ = new T_JacobianForm(V, V);
        T_ResidualForm* residualForm_ = new T_ResidualForm(V);
        _a->u = *_u;
        _a->lmbda = lambda; _a->dt = dt; _a->theta = theta;
        _L->u = *_u; _L->u0 = *_u0;
        _L->lmbda = lambda; _L->dt = dt; _L->theta = theta;
  
        // Wrap pointers in a smart pointer
        a.reset(_a);
        L.reset(_L);
  
        // Set solution to intitial condition
        InitialConditions u_initial;
        *_u = u_initial;
  */
      }

      virtual void processAdditionalVector (dolfin::Vector& vec);
  
      // Function space, forms and functions
      std::shared_ptr<dolfin::FunctionSpace> functionSpace_;
      std::shared_ptr<dolfin::FunctionSpace> additionalFunctionSpace_;
      T_ResidualForm residualForm_;
      T_JacobianForm jacobianForm_;
      T_AdditionalForm additionalForm_;
      std::shared_ptr<dolfin::Function> solution_;
      std::shared_ptr<dolfin::Function> previousSolution_;
      std::map <std::string, dolfin::DirichletBC> dirichletBCs_;
  };



    // ==============================================================================================//
    // ==================================== IMPLEMENTATION ==========================================//
    // ==============================================================================================//



  template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_ResidualForm, class T_JacobianForm, class T_AdditionalForm>
        void UflToNewton <T_FunctionSpace, T_AdditionalFunctionSpace, T_ResidualForm, T_JacobianForm, T_AdditionalForm>::
        setCoefficient (const std::string& coefficientType, 
                        const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                        const std::string& coefficientName)
        {
            if (coefficientType == "residual_form")
            {
                dolfin::log (dolfin::DBG, "Setting residual form coefficient \"%s\"...", coefficientName.c_str ());
                residualForm_.set_coefficient (coefficientName, coefficientValue);
            }
            else if (coefficientType == "jacobian_form")
            {
                dolfin::log (dolfin::DBG, "Setting jacobian form coefficient \"%s\"...", coefficientName.c_str ());
                jacobianForm_.set_coefficient (coefficientName, coefficientValue);
            }
/*            else if (coefficientType == "initial_guess")
            {
                dolfin::log (dolfin::DBG, "Setting initial guess...");
                // check whether coefficientValue is a pointer to dolfin::Function or dolfin::Expression
                if (std::dynamic_pointer_cast<const dolfin::Function> (coefficientValue) != nullptr)
                {
                    solution_.back ().second = *(std::dynamic_pointer_cast<const dolfin::Function> (coefficientValue));
                }
                else if (std::dynamic_pointer_cast<const dolfin::Expression> (coefficientValue) != nullptr)
                {
                    solution_.back ().second = *(std::dynamic_pointer_cast<const dolfin::Expression> (coefficientValue));
                }
                else
                {
                    dolfin::warning ("Cannot set initial guess in nonlinear differential problem. Input argument is neither a dolfin::Function nor a dolfin::Expression");
                }
            }*/
            else if (coefficientType == "additional_form")
            {
                dolfin::log (dolfin::DBG, "Setting jacobian form coefficient \"%s\"...", coefficientName.c_str ());
                additionalForm_.set_coefficient (coefficientName, coefficientValue);
            }
            else
            {
                dolfin::warning ("Cannot set coefficient in non linear differential problem. Coefficient type \"%s\" unknown",
                                 coefficientType.c_str ());
            }
        }

  template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_ResidualForm, class T_JacobianForm, class T_AdditionalForm>
        void UflToNewton <T_FunctionSpace, T_AdditionalFunctionSpace, T_ResidualForm, T_JacobianForm, T_AdditionalForm>::
        setIntegrationSubdomain (const std::string& formType,
                                  std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                  const dcp::SubdomainType& subdomainType)
        {
            if (formType == "residual_form")
            {
                if (subdomainType == dcp::SubdomainType::INTERNAL_CELLS)
                {
                    dolfin::log (dolfin::DBG, "Setting residual form integration subdomain on INTERNAL_CELLS...");
                    residualForm_.set_cell_domains (meshFunction);
                }
                else if (subdomainType == dcp::SubdomainType::INTERNAL_FACETS)
                {
                    dolfin::log (dolfin::DBG, "Setting residual form integration subdomain on INTERNAL_FACETS...");
                    residualForm_.set_interior_facet_domains (meshFunction);
                }
                else if (subdomainType == dcp::SubdomainType::BOUNDARY_FACETS)
                {
                    dolfin::log (dolfin::DBG, "Setting residual form integration subdomain on BOUNDARY_FACETS...");
                    residualForm_.set_exterior_facet_domains (meshFunction);
                }
                else if (subdomainType == dcp::SubdomainType::VERTICES)
                {
                    dolfin::log (dolfin::DBG, "Setting residual form integration subdomain on VERTICES...");
                    residualForm_.set_vertex_domains (meshFunction);
                }
                else
                {
                    dolfin::warning ("unknown subdomain type requested while trying to apply mesh function to residual form"); 
                }
            }
            else if (formType == "jacobian_form")
            {
                if (subdomainType == dcp::SubdomainType::INTERNAL_CELLS)
                {
                    dolfin::log (dolfin::DBG, "Setting jacobian form integration subdomain on INTERNAL_CELLS...");
                    jacobianForm_.set_cell_domains (meshFunction);
                }
                else if (subdomainType == dcp::SubdomainType::INTERNAL_FACETS)
                {
                    dolfin::log (dolfin::DBG, "Setting jacobian form integration subdomain on INTERNAL_FACETS...");
                    jacobianForm_.set_interior_facet_domains (meshFunction);
                }
                else if (subdomainType == dcp::SubdomainType::BOUNDARY_FACETS)
                {
                    dolfin::log (dolfin::DBG, "Setting jacobian form integration subdomain on BOUNDARY_FACETS...");
                    jacobianForm_.set_exterior_facet_domains (meshFunction);
                }
                else if (subdomainType == dcp::SubdomainType::VERTICES)
                {
                    dolfin::log (dolfin::DBG, "Setting jacobian form integration subdomain on VERTICES...");
                    jacobianForm_.set_vertex_domains (meshFunction);
                }
                else
                {
                    dolfin::warning ("unknown subdomain type requested while trying to apply mesh function to jacobian form"); 
                }
            }
            else if (formType == "additional_form")
            {
                if (subdomainType == dcp::SubdomainType::INTERNAL_CELLS)
                {
                    dolfin::log (dolfin::DBG, "Setting jacobian form integration subdomain on INTERNAL_CELLS...");
                    additionalForm_.set_cell_domains (meshFunction);
                }
                else if (subdomainType == dcp::SubdomainType::INTERNAL_FACETS)
                {
                    dolfin::log (dolfin::DBG, "Setting jacobian form integration subdomain on INTERNAL_FACETS...");
                    additionalForm_.set_interior_facet_domains (meshFunction);
                }
                else if (subdomainType == dcp::SubdomainType::BOUNDARY_FACETS)
                {
                    dolfin::log (dolfin::DBG, "Setting jacobian form integration subdomain on BOUNDARY_FACETS...");
                    additionalForm_.set_exterior_facet_domains (meshFunction);
                }
                else if (subdomainType == dcp::SubdomainType::VERTICES)
                {
                    dolfin::log (dolfin::DBG, "Setting jacobian form integration subdomain on VERTICES...");
                    additionalForm_.set_vertex_domains (meshFunction);
                }
                else
                {
                    dolfin::warning ("unknown subdomain type requested while trying to apply mesh function to jacobian form"); 
                }
            }
            else
            {
                dolfin::warning ("Cannot set integration subdomain in linear differential problem. Form type \"%s\" unknown",
                                 formType.c_str ());
            }

        }

  template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_ResidualForm, class T_JacobianForm, class T_AdditionalForm>
        void UflToNewton <T_FunctionSpace, T_AdditionalFunctionSpace, T_ResidualForm, T_JacobianForm, T_AdditionalForm>::
        processAdditionalVector (dolfin::Vector& vec)
        {
            //const dolfin::MeshFunction<std::size_t>& meshFunction (* additionalForm_.vertex_domains());
            auto& meshFunctionVertex (* additionalForm_.vertex_domains());
            auto& meshFunctionFacet (* additionalForm_.exterior_facet_domains());
/*            const dolfin::Mesh& mesh (* meshFunction.mesh());
std::cerr << ".h OK fino a " << __LINE__ << std::endl;
!!! GIVES SEGMENTATION FAULT*/
            const dolfin::Mesh& mesh (* additionalFunctionSpace_->mesh());
            additionalFunctionSpace_->print_dofmap();
(*additionalFunctionSpace_)[0]->print_dofmap();
(*(*additionalFunctionSpace_)[0])[0]->print_dofmap();
(*additionalFunctionSpace_)[1]->print_dofmap();
std::cerr << ".h OK fino a " << __LINE__ << std::endl;
            const dolfin::GenericDofMap& dofMap (* additionalFunctionSpace_->dofmap());
std::cerr << "dofmap dimension " << dofMap.geometric_dimension() << std::endl;
            std::vector<double> coordinates (dofMap.tabulate_all_coordinates(mesh));
std::cerr << "coordinates size " << coordinates.size() << std::endl;
std::ostream_iterator<double> out_it (std::cerr,", ");
std::ostream_iterator<std::size_t> out_it_int (std::cerr,", ");
std::copy ( coordinates.begin(), coordinates.end(), out_it );
            std::vector<dolfin::la_index> dofs (dofMap.dofs());
std::cerr << std::endl << "dofs size " << dofs.size() << std::endl;
std::copy ( dofs.begin(), dofs.end(), out_it_int );
/* !!! SEGFAULT IN WHAT FOLLOWS
            const dolfin::GenericDofMap& velocityDofMap (* dofMap.extract_sub_dofmap({0},mesh));
            const dolfin::GenericDofMap& velocityDofMap (* (*additionalFunctionSpace_)[0]->dofmap());
std::cerr << "\nVELOCITY" << std::endl;
            coordinates = velocityDofMap.tabulate_all_coordinates(mesh);
std::cerr << "coordinates size " << coordinates.size() << std::endl;
std::copy ( coordinates.begin(), coordinates.end(), out_it );
            dofs = velocityDofMap.dofs();
std::cerr << std::endl << "dofs size " << dofs.size() << std::endl;
            const dolfin::GenericDofMap& pressureDofMap (* dofMap.extract_sub_dofmap({1},mesh));
std::cerr << "\nPRESSURE" << std::endl;
            coordinates = pressureDofMap.tabulate_all_coordinates(mesh);
std::cerr << "coordinates size " << coordinates.size() << std::endl;
std::copy ( coordinates.begin(), coordinates.end(), out_it );
            dofs = pressureDofMap.dofs();
std::cerr << std::endl << "dofs size " << dofs.size() << std::endl;*/
dolfin::Function xy (* additionalFunctionSpace_);
xy = XY();
const dolfin::GenericVector& xyvec (* xy.vector());
std::cerr << std::endl << "XY function vector (size = " << xyvec.size() << ")" << std::endl;
for (std::size_t i=0; i!= xyvec.size(); ++i) std::cerr << xyvec[i] << ", ";
const std::vector<double>& xx (mesh.geometry().x());
std::cerr << std::endl << "xx vector (size = " << xx.size() << ")" << std::endl;
std::copy ( xx.begin(), xx.end(), out_it );
if (! meshFunctionVertex.empty())
  {std::cerr << std::endl << "mesh function vertex values (size = " << meshFunctionVertex.size() << ")" << std::endl;
  for (std::size_t i=0; i!= meshFunctionVertex.size(); ++i) std::cerr << meshFunctionVertex[i] << ", ";}
else std::cerr << std::endl << "meshFunctionVertex is empty" << std::endl;
if (! meshFunctionFacet.empty())
  {std::cerr << std::endl << "mesh function facet values (size = " << meshFunctionFacet.size() << ")" << std::endl;
  for (std::size_t i=0; i!= meshFunctionFacet.size(); ++i) std::cerr << meshFunctionFacet[i] << ", ";}
else std::cerr << std::endl << "meshFunctionFacet is empty" << std::endl;
std::cerr << std::endl 
          << std::endl
          << std::endl;
const dolfin::GenericDofMap& dofmapComponent (*(*additionalFunctionSpace_)[0]->dofmap());
const std::vector<unsigned int>& meshcells (mesh.cells());
std::cerr << "mesh cells" << std::endl; std::copy (meshcells.begin(), meshcells.end(), out_it); std::cerr << std::endl;
/*std::vector<double> dmCompCoords (dofmapComponent.tabulate_all_coordinates(mesh));
std::cerr << "component coordinates" << std::endl; std::copy (dmCompCoords.begin(), dmCompCoords.end(), out_it); std::cerr << std::endl;*/
            std::vector<dolfin::la_index> dmCompDofs (dofmapComponent.dofs());
std::cerr << "component dofs (size " << dmCompDofs.size() << ")" << std::endl; std::copy (dmCompDofs.begin(), dmCompDofs.end(), out_it); std::cerr << std::endl;
std::vector<dolfin::la_index> okDofs;
dolfin::Array<double> pointCoords (2);
TriplePointLeftVertex tplv;
TriplePointRightVertex tprv;
for (auto it=dmCompDofs.begin(); it!=dmCompDofs.end(); it++)
{
  pointCoords[0] = coordinates[2*(*it)];
  pointCoords[1] = coordinates[2*(*it)+1];
  if (! (tplv.inside(pointCoords,true) || tprv.inside(pointCoords,true)))
  { 
    okDofs.push_back(*it);
    const double val[1] = {0};
    const dolfin::la_index idx[1] = {*it};
    vec.set(val,1,idx);
  }
  else
    std::cerr << pointCoords[0] << "  " << pointCoords[1] << std::endl;
}
std::cerr << "selected dofs (size " << okDofs.size() << ")" << std::endl;
std::copy (okDofs.begin(), okDofs.end(), out_it_int);
std::cerr << std::endl;

/*for (auto it = meshcells.begin(); it!=meshcells.end(); it++)
{
  
}
for (std::size_t cellIdx = 0; cellIdx!=mesh.size(dofmapComponent.geometric_dimension()); ++cellIdx)
{
  const std::vector<dolfin::la_index>& cellDofs (dofmapComponent.cell_dofs(cellIdx));
std::cerr << "cell " << cellIdx << std::endl;
std::copy (cellDofs.begin(),cellDofs.end(),out_it);
std::cerr << std::endl;
  boost::multi_array<double,2> coords;
  double vCoords[3*2];
  dolfin::Cell cell (mesh, cellIdx);
  cell.get_vertex_coordinates (vCoords);
  std::vector<double> vertexCoords (vCoords, vCoords + sizeof(vCoords)/sizeof(double));
  dofmapComponent.tabulate_coordinates (coords,vertexCoords,cell);
std::cerr << "coords " << std::endl;
std::copy (coords.begin(),coords.end(),out_it);
std::cerr << std::endl;
std::cerr << "vertexCoords " << std::endl;
//for (std::size_t i=0; i!=2; ++i) for (std::size_t j=0; j!=coords.size()/2; ++j) std::cerr << coords[i][j] << ", ";
std::cerr << coords[0][1] << " ";
std::cerr << std::endl;
}*/
        }

} //end of namespace
#endif
