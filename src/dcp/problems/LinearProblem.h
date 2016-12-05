/* 
 *  Copyright (C) 2014, Mattia Tamellini, mattia.tamellini@gmail.com
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

#ifndef SRC_PROBLEMS_LINEARPROBLEM_H_INCLUDE_GUARD
#define SRC_PROBLEMS_LINEARPROBLEM_H_INCLUDE_GUARD

#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/la/LinearSolver.h>
#include <dolfin/la/GenericLinearSolver.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/parameter/Parameters.h>
#include <vector>
#include <string>
#include <memory>
#include <dcp/problems/GenericLinearProblem.h>
#include <dcp/factories/LinearSolverFactory.h>
#include <dcp/problems/SubdomainType.h>

namespace dcp
{
    /*! \class LinearProblem LinearProblem.h
     *  \brief Class for linear differential problems.
     *
     *  This class represents problem of the form
     *  \f[
     *      \mbox{Find } u \in V : a \left(u, v\right) = F \left(v\right) \ \forall\,v\,\in\,V
     *  \f]
     *  with \f$ a \left(u, v\right) : V \times V \rightarrow \mathds{R}\f$ bilinear form on \f$V\f$
     *  and \f$ L \left(v\right) : V \rightarrow \mathds{R} \f$ linear form on the same space.
     *  
     *  It inherits publicly from \c GenericLinearProblem and it extends its functionalities to a concrete 
     *  differential problem.
     *  Template arguments are:
     *  \arg T_BilinearForm the bilinear form type
     *  \arg T_LinearForm the linear form type
     *  \arg T_LinearSolverFactory the type of the factory that creates the linear solver. By default, it is set
     *  to \c dcp::LinearSolverFactory
     */

    template <class T_BilinearForm_, class T_LinearForm_, class T_LinearSolverFactory_ = dcp::LinearSolverFactory>
        class LinearProblem : public dcp::GenericLinearProblem
        {
            // ---------------------------------------------------------------------------------------------//  

            public:
                /******************* TYPEDEFS *******************/
                typedef T_BilinearForm_        T_BilinearForm;
                typedef T_LinearForm_          T_LinearForm;
                typedef T_LinearSolverFactory_ T_LinearSolverFactory;
                
                
                /******************* CONSTRUCTORS *******************/
                //! Default constructor is deleted. The class is not default constructable.
                LinearProblem () = delete;

                //!  Constructor [1]
                /*!
                 *  \param functionSpace the problem finite element space 
                 *  The constructors also sets the following parameters:
                 *      - \c "lump_matrix" boolean value that controls whether the system matrix is lumped or not.
                 *        Default value: false
                 *      - \c "solver_type" the type of the solver to be used. Default value: \c "lu_solver"
                 *      - \c "solver_method" the method to be used by the solver. Default value: \c "default"
                 *      - \c "solver_preconditioner" the preconditioner to be used by the solver. 
                 *        Default value: \c "default"
                 *      - \c "system_is_assembled" a flag that can be set to false if one wants to force a reassembly of
                 *        the linear system when the \ solve method is called. Default value: \c false
                 *      - \c "force_reassemble_system" a flag that, if set to \c true, causes the system to be 
                 *        reassembled every time the \c solve method is called. Default value: \c false
                 */
                LinearProblem (const std::shared_ptr<const dolfin::FunctionSpace> functionSpace);

                //!  Constructor [2]
                /*!
                 *  \param functionSpace the problem finite element space as a <tt>const std::shared_ptr </tt> to 
                 *  \c dolfin::FunctionSpace
                 *  \param bilinearForm a \c const reference to the problem's bilinear form
                 *  \param linearForm a \c const reference to the problem's linear form
                 *  The stored function space's ownership will be shared between the object and the input argument.  The
                 *  bilinear and linear form will be created too, calling the constructor which takes the function space
                 *  as input.
                 *  The constructors also sets the following parameters:
                 *      - \c "lump_matrix" boolean value that controls whether the system matrix is lumped or not.
                 *        Default value: false
                 *      - \c "solver_type" the type of the solver to be used. Default value: \c "lu_solver"
                 *      - \c "solver_method" the method to be used by the solver. Default value: \c "default"
                 *      - \c "solver_preconditioner" the preconditioner to be used by the solver. 
                 *        Default value: \c "default"
                 *      - \c "system_is_assembled" a flag that can be set to false if one wants to force a reassembly of
                 *        the linear system when the \ solve method is called. Default value: \c false
                 *      - \c "force_reassemble_system" a flag that, if set to \c true, causes the system to be 
                 *        reassembled every time the \c solve method is called. Default value: \c false
                 */
                LinearProblem (const std::shared_ptr<const dolfin::FunctionSpace> functionSpace,
                               const T_BilinearForm& bilinearForm,
                               const T_LinearForm& linearForm);


                /******************* DESTRUCTOR *******************/

                //! Destructor
                /*! 
                 *  Default destructor, since members of the class are trivially 
                 *  destructible.
                 */
                virtual ~LinearProblem () {};

                
                /******************* GETTERS *******************/
                //! Get const reference to the problem's linear form
                /*! 
                 *  \return a const reference to the problem's linear form
                 */
                virtual const T_BilinearForm& bilinearForm () const override;

                //! Get const reference to the problem's linear form
                /*! 
                 *  \return a const reference to the problem's linear form
                 */
                virtual const T_LinearForm& linearForm () const override;

                //! Get const reference to the problem's linear operator
                /*!
                 *  \return a const reference to the problem's linear operator, which
                 *  is a \c dolfin::Matrix
                 */
                virtual const dolfin::Matrix& linearOperator () const override;

                //! Get const reference to the problem's right hand side
                /*!
                 *  \return a const reference to the problem's right hand side, which
                 *  is a \c dolfin::Vector
                 */
                virtual const dolfin::Vector& rhs () const override;


                /******************* SETTERS *******************/
                
                //! Set coefficient [1]. Override of virtual function in \c GenericProblem.
                /*!
                 *  Possible values for \c coefficientType are:
                 *  \li \c bilinear_form to set the coefficient in the bilinear form
                 *  \li \c linear_form to set the coefficient in the linear form
                 *  
                 *  See \c GenericProblem documentation for more details on the function
                 */
                virtual void setCoefficient (const std::string& coefficientType, 
                                             const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                             const std::string& coefficientName) override;

                //! Set coefficient [2]. Override of virtual function in \c GenericProblem.
                /*!
                 *  Possible values for \c coefficientType are:
                 *  \li \c bilinear_form to set the coefficient in the bilinear form
                 *  \li \c linear_form to set the coefficient in the linear form
                 *  
                 *  See \c GenericProblem documentation for more details on the function
                 */
                virtual void setCoefficient (const std::string& coefficientType,
                                             const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                             const std::size_t& coefficientNumber) override;

                //! Set integration subdomains for the forms. Override of virtual function in \c GenericProblem
                /*! 
                 *  Possible values for \c formType are:
                 *  \li \c bilinear_form to set the integration subdomain in the bilinear form
                 *  \li \c linear_form to set the integration subdomain in the linear form
                 *  
                 *  See \c GenericProblem documentation for more details on the function
                 */
                virtual void setIntegrationSubdomain (const std::string& formType,
                                                      std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                                      const dcp::SubdomainType& subdomainType) override;

                //! Add Dirichlet boundary condition to the problem [1]. Overrides method in \c GenericProblem
                /*!
                 *  This method adds to the base class method the setting of parameter \c system_is_assembled to 
                 *  \c false.
                 *  \param condition the boundary condition to enforce
                 *  \param boundary the boundary on which to enforce the condition
                 *  \param bcName the name identifying the boundary condition. If empty, 
                 *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
                 *  
                 *  \return boolean flag, with \c true representing success and \c false representing failure
                 */
                virtual bool addDirichletBC (const dolfin::GenericFunction& condition, 
                                             const dolfin::SubDomain& boundary,
                                             std::string bcName = "") override;

                //! Add Dirichlet boundary condition to the problem [2]. Overrides method in \c GenericProblem
                /*!
                 *  This method adds to the base class method the setting of parameter \c system_is_assembled to 
                 *  \c false.
                 *  \param condition the boundary condition to enforce
                 *  \param boundary the boundary on which to enforce the condition
                 *  \param component the function space component on which the boundary condition should be imposed. 
                 *  \param bcName the name identifying the boundary condition. If empty, 
                 *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
                 *  
                 *  \return boolean flag, with \c true representing success and \c false representing failure
                 */
                virtual bool addDirichletBC (const dolfin::GenericFunction& condition, 
                                             const dolfin::SubDomain& boundary,
                                             const std::size_t& component,
                                             std::string bcName = "") override;

                //! Add Dirichlet boundary condition to the problem [3]
                /*!
                 *  This method adds to the base class method the setting of parameter \c system_is_assembled to 
                 *  \c false.
                 *  \param condition the boundary condition to enforce
                 *  \param boundary the boundary on which to enforce the condition
                 *  \param bcName the name identifying the boundary condition. If empty,
                 *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
                 *  
                 *  \return boolean flag, with \c true representing success and \c false representing failure
                 */
                virtual bool addDirichletBC (std::shared_ptr<const dolfin::GenericFunction> condition, 
                                             std::shared_ptr<const dolfin::SubDomain> boundary,
                                             std::string bcName = "") override;

                //! Add Dirichlet boundary condition to the problem [4]
                /*!
                 *  This method adds to the base class method the setting of parameter \c system_is_assembled to 
                 *  \c false.
                 *  \param condition the boundary condition to enforce
                 *  \param boundary the boundary on which to enforce the condition
                 *  \param component the function space component on which the boundary condition should be imposed. 
                 *  For instance, this can be useful if we have a vector space and we want only the orizontal component to
                 *  have a fixed value
                 *  \param bcName the name identifying the boundary condition. If empty,
                 *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
                 *  
                 *  \return boolean flag, with \c true representing success and \c false representing failure
                 */
                virtual bool addDirichletBC (std::shared_ptr<const dolfin::GenericFunction> condition, 
                                             std::shared_ptr<const dolfin::SubDomain> boundary,
                                             const std::size_t& component,
                                             std::string bcName = "") override;

                //! Add Dirichlet boundary condition to the problem [5]. Overrides method in \c GenericProblem
                /*!
                 *  This method adds to the base class method the setting of parameter \c system_is_assembled to \c false.
                 *  \param dirichletCondition a const reference to the dirichlet boundary condition to be added to the
                 *  problem
                 *  \param bcName the name identifying the boundary condition. If empty, 
                 *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
                 *  
                 *  \return boolean flag, with \c true representing success and \c false representing failure
                 */
                virtual bool addDirichletBC (const dolfin::DirichletBC& dirichletCondition, 
                                             std::string bcName = "") override;

                //! Add Dirichlet boundary condition to the problem [6]. Overrides method in \c GenericProblem
                /*!
                 *  This method adds to the base class method the setting of parameter \c system_is_assembled to \c false.
                 *  \param dirichletCondition a rvalue reference to the dirichlet boundary condition to be added to the
                 *  problem
                 *  \param bcName the name identifying the boundary condition. If empty, 
                 *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
                 *  
                 *  \return boolean flag, with \c true representing success and \c false representing failure
                 */
                virtual bool addDirichletBC (dolfin::DirichletBC&& dirichletCondition, 
                                             std::string bcName = "") override;

                //! Remove Dirichlet boundary condition with given position. Overrides method in \c GenericProblem
                /*!
                 *  This method adds to the base class method the setting of parameter \c system_is_assembled to 
                 *  \c false.
                 *  \param bcName name of the boundary condition to be removed.
                 *  
                 *  \return boolean flag, with \c true representing success and \c false representing failure
                 */
                virtual bool removeDirichletBC (const std::string& bcName) override;
                
                //! Method to update class members. It checks for differences between desired and current solver
                //! parameters and creates a new solver, setting also the proper parameters
                virtual void update () override;

                
                /******************* METHODS *******************/
                //! Lump system matrix
                /*! 
                 *  Performs lumping of the system matrix by substituting the diagonal with the sum of each row and
                 *  setting the extra-diagonal terms to zero.
                 */
                virtual void lumpMatrix () override;

                //! Assemble the linear system
                /*!
                 *  Assemble the linear system related to the problem to be solved
                 */
                virtual void assembleLinearSystem () override;

                //! Solve problem
                /*!
                 *  This method solves the problem defined. It uses the private members' value to set the problem and then
                 *  stores the solution in the private member \c solution_. See documentation of \c dcp::GenericProblem
                 *  for more details on how the protected member \c solution_ works and why it is declared as a 
                 *  \c std::pair.
                 *  Note that the method \c solve() checks the problem's member
                 *  \c parameters to decided whether problem's matrix and vector should be reassembled and if
                 *  the values of the parameters \c solver_type, \c solver_method and \c solver_preconditioner 
                 *  match the values of \c solverType_, \c solverMethod_ and \c solverPreconditioner_.
                 *  If they differ, it calls \c createSolver_().
                 * 
                 *  \param solveType the solution type wanted. Possible values are:
                 *  \li \c "default" : the normal solution process
                 *  \li \c "stash" : solve the problem but store the solution in \c stashedSolution_
                 */
                virtual void solve (const std::string& solveType = "default") override;

                //! Clone method. Overrides method in \c GenericProblem
                /*!
                 *  It uses the parameter \c clone_method to decide which type of cloning to perform.
                 *  Possible values for such parameter are:
                 *  \li deep_clone the new object is created calling the constructor that takes a mesh and a function 
                 *  space as input, thus creating a copy of such objects and returning a completely independent 
                 *  cloned object. 
                 *  \li shallow_clone calls the constructor that takes shared pointers as input: the mesh and
                 *  the function space are not copied but shared between the current object and its clone. 
                 *  
                 *  The default value for parameter \c clone_method is \c shallow_clone
                 *  
                 *  \return a pointer to the cloned object
                 */
                virtual dcp::LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>* clone () const override;

                // ---------------------------------------------------------------------------------------------//

            protected:
                //! Creates a linear solver object of type passed as input. 
                /*!
                 *  The solver will be created using the parameters \c solver_type, \c solver_method
                 *  and \c solver_preconditioner. It will also set the protected members \c solverType_, \c
                 *  solverMethod_ and \c solverPreconditioner_.
                 */
                std::unique_ptr<dolfin::GenericLinearSolver> createSolver_ ();
                
                //! The bilinear form
                T_BilinearForm bilinearForm_;

                //! The linear form
                T_LinearForm linearForm_;

                //! The solver
                /*! 
                 *  We use a pointer so that polymorphism can be applied.
                 */
                std::unique_ptr<dolfin::GenericLinearSolver> solver_;
                
                //! Matrix to hold the problem's discrete operator
                dolfin::Matrix problemMatrix_;

                //! Vector to store the right hand side of the discrete problem
                dolfin::Vector rhsVector_;

                //! String to contain the current solver type being used
                std::string solverType_;

                //! String to contain the current solver method being used
                std::string solverMethod_;

                //! String to contain the current solver preconditioner being used
                std::string solverPreconditioner_;
                // ---------------------------------------------------------------------------------------------//

            private:
        };



    // ============================================================================================== //
    // ==================================== IMPLEMENTATION ========================================== //
    // ============================================================================================== //

    /******************* CONSTRUCTORS *******************/

    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        LinearProblem (const std::shared_ptr<const dolfin::FunctionSpace> functionSpace) :
            GenericLinearProblem (functionSpace),
            bilinearForm_ (functionSpace_, functionSpace_),
            linearForm_ (functionSpace_),
            solver_ (nullptr),
            problemMatrix_ (),
            rhsVector_ (),
            solverType_ ("lu_solver"),
            solverMethod_ ("default"),
            solverPreconditioner_ ("default")
        { 
            dolfin::begin (dolfin::DBG, "Building LinearProblem...");

            solution_.emplace_back (std::make_pair (-1, dolfin::Function (functionSpace_)));

            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("lump_matrix", false);
            parameters.add ("solver_type", solverType_);
            parameters.add ("solver_method", solverMethod_);
            parameters.add ("solver_preconditioner", solverPreconditioner_);
            parameters.add ("system_is_assembled", false);
            parameters.add ("force_reassemble_system", false);  

            dolfin::begin (dolfin::DBG, "Creating solver...");
            solver_ = createSolver_ ();
            dolfin::end ();

            dolfin::end ();

            dolfin::log (dolfin::DBG, "LinearProblem object created");
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        LinearProblem (const std::shared_ptr<const dolfin::FunctionSpace> functionSpace,
                       const T_BilinearForm& bilinearForm,
                       const T_LinearForm& linearForm) :
            GenericLinearProblem (functionSpace),
            bilinearForm_ (bilinearForm),
            linearForm_ (linearForm),
            solver_ (nullptr),
            problemMatrix_ (),
            rhsVector_ (),
            solverType_ ("lu_solver"),
            solverMethod_ ("default"),
            solverPreconditioner_ ("default")
        { 
            dolfin::begin (dolfin::DBG, "Building LinearProblem...");
            
            solution_.emplace_back (std::make_pair (-1, dolfin::Function (functionSpace_)));
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("lump_matrix", false);
            parameters.add ("solver_type", solverType_);
            parameters.add ("solver_method", solverMethod_);
            parameters.add ("solver_preconditioner", solverPreconditioner_);
            parameters.add ("system_is_assembled", false);
            parameters.add ("force_reassemble_system", false);
            
            dolfin::begin (dolfin::DBG, "Creating solver...");
            solver_ = createSolver_ ();
            dolfin::end ();
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "LinearProblem object created");
        }

    

    /******************* GETTERS *******************/

    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        const T_BilinearForm& LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        bilinearForm () const
        {
            return bilinearForm_;
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        const T_LinearForm& LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        linearForm () const
        {
            return linearForm_;
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        const dolfin::Matrix& LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        linearOperator () const
        {
            return problemMatrix_;
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        const dolfin::Vector& LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        rhs () const
        {
            return rhsVector_;
        }



    /******************* SETTERS *******************/

    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        void LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        setCoefficient (const std::string& coefficientType, 
                        const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                        const std::string& coefficientName)
        {
            if (coefficientType == "bilinear_form")
            {
                dolfin::log (dolfin::DBG, "Setting bilinear form coefficient \"%s\"...", coefficientName.c_str ());
                bilinearForm_.set_coefficient (coefficientName, coefficientValue);
                parameters ["system_is_assembled"] = false;
            }
            else if (coefficientType == "linear_form")
            {
                dolfin::log (dolfin::DBG, "Setting linear form coefficient \"%s\"...", coefficientName.c_str ());
                linearForm_.set_coefficient (coefficientName, coefficientValue);
                parameters ["system_is_assembled"] = false;
            }
            else
            {
                dolfin::warning
                    ("Cannot set coefficient \"%s\" in linear differential problem. Form type \"%s\" unknown",
                     coefficientName.c_str (),
                     coefficientType.c_str ());
            }

        }

    

    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        void LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        setCoefficient (const std::string& coefficientType,
                        const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                        const std::size_t& coefficientNumber)
        {
            if (coefficientType == "bilinear_form")
            {
                dolfin::log (dolfin::DBG, "Setting bilinear form coefficient number %d...", coefficientNumber);
                bilinearForm_.set_coefficient (coefficientNumber, coefficientValue);
                parameters ["system_is_assembled"] = false;
            }
            else if (coefficientType == "linear_form")
            {
                dolfin::log (dolfin::DBG, "Setting linear form coefficient number %d...", coefficientNumber);
                linearForm_.set_coefficient (coefficientNumber, coefficientValue);
                parameters ["system_is_assembled"] = false;
            }
            else
            {
                dolfin::warning
                    ("Cannot set coefficient number %d in linear differential problem. Form type \"%s\" unknown",
                     coefficientNumber,
                     coefficientType.c_str ());
            }
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        void LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        setIntegrationSubdomain (const std::string& formType,
                                  std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                  const dcp::SubdomainType& subdomainType)
        {
            if (formType == "bilinear_form")
            {
                if (subdomainType == dcp::SubdomainType::INTERNAL_CELLS)
                {
                    dolfin::log (dolfin::DBG, "Setting bilinear form integration subdomain on INTERNAL_CELLS...");
                    bilinearForm_.set_cell_domains (meshFunction);
                    parameters ["system_is_assembled"] = false;
                }
                else if (subdomainType == dcp::SubdomainType::INTERNAL_FACETS)
                {
                    dolfin::log (dolfin::DBG, "Setting bilinear form integration subdomain on INTERNAL_FACETS...");
                    bilinearForm_.set_interior_facet_domains (meshFunction);
                    parameters ["system_is_assembled"] = false;
                }
                else if (subdomainType == dcp::SubdomainType::BOUNDARY_FACETS)
                {
                    dolfin::log (dolfin::DBG, "Setting bilinear form integration subdomain on BOUNDARY_FACETS...");
                    bilinearForm_.set_exterior_facet_domains (meshFunction);
                    parameters ["system_is_assembled"] = false;
                }
                else
                {
                    dolfin::warning ("unknown subdomain type requested while trying to apply mesh function to bilinear form"); 
                }
            }
            else if (formType == "linear_form")
            {
                if (subdomainType == dcp::SubdomainType::INTERNAL_CELLS)
                {
                    dolfin::log (dolfin::DBG, "Setting linear form integration subdomain on INTERNAL_CELLS...");
                    linearForm_.set_cell_domains (meshFunction);
                    parameters ["system_is_assembled"] = false;
                }
                else if (subdomainType == dcp::SubdomainType::INTERNAL_FACETS)
                {
                    dolfin::log (dolfin::DBG, "Setting linear form integration subdomain on INTERNAL_FACETS...");
                    linearForm_.set_interior_facet_domains (meshFunction);
                    parameters ["system_is_assembled"] = false;
                }
                else if (subdomainType == dcp::SubdomainType::BOUNDARY_FACETS)
                {
                    dolfin::log (dolfin::DBG, "Setting linear form integration subdomain on BOUNDARY_FACETS...");
                    linearForm_.set_exterior_facet_domains (meshFunction);
                }
                else
                {
                    dolfin::warning ("unknown subdomain type requested while trying to apply mesh function to linear form"); 
                }
            }
            else
            {
                dolfin::warning ("Cannot set integration subdomain in linear differential problem. Form type \"%s\" unknown",
                                 formType.c_str ());
            }

        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        bool LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        addDirichletBC (const dolfin::GenericFunction& condition, 
                        const dolfin::SubDomain& boundary,
                        std::string bcName)
        {
            return addDirichletBC (dolfin::reference_to_no_delete_pointer (condition), 
                                   dolfin::reference_to_no_delete_pointer (boundary),
                                   bcName); 
        }
    


    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        bool LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        addDirichletBC (const dolfin::GenericFunction& condition, 
                        const dolfin::SubDomain& boundary,
                        const std::size_t& component,
                        std::string bcName)
        {
            return addDirichletBC (dolfin::reference_to_no_delete_pointer (condition), 
                                   dolfin::reference_to_no_delete_pointer (boundary),
                                   component,
                                   bcName); 
        }
                    


    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        bool LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        addDirichletBC (std::shared_ptr<const dolfin::GenericFunction> condition, 
                        std::shared_ptr<const dolfin::SubDomain> boundary,
                        std::string bcName)
        {
            if (bcName.empty ())
            {
                bcName = "dirichlet_condition_" + std::to_string (dirichletBCsCounter_);
                dirichletBCsCounter_++;
            }

            dolfin::log (dolfin::DBG, 
                         "Adding dirichlet boundary condition to boundary conditions map with name \"%s\"...",
                         bcName.c_str ());

            auto result = dirichletBCs_.emplace (bcName, dolfin::DirichletBC (functionSpace_, condition, boundary));

            if (result.second == false)
            {
                dolfin::warning ("DirichletBC object not inserted because key \"%s\" already in map",
                                 bcName.c_str ());
            }
            else
            {
                parameters ["system_is_assembled"] = false;
            }

            return result.second;
        }
    


    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        bool LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        addDirichletBC (std::shared_ptr<const dolfin::GenericFunction> condition, 
                        std::shared_ptr<const dolfin::SubDomain> boundary,
                        const std::size_t& component,
                        std::string bcName)
        {
            if (bcName.empty ())
            {
                bcName = "dirichlet_condition_" + std::to_string (dirichletBCsCounter_);
                dirichletBCsCounter_++;
            }

            dolfin::log (dolfin::DBG, 
                         "Adding dirichlet boundary condition to boundary conditions map with name \"%s\"...",
                         bcName.c_str ());

            auto result = dirichletBCs_.emplace
                (bcName, dolfin::DirichletBC ((*functionSpace_) [component], condition, boundary));

            if (result.second == false)
            {
                dolfin::warning ("DirichletBC object not inserted because key \"%s\" already in map",
                                 bcName.c_str ());
            }
            else
            {
                parameters ["system_is_assembled"] = false;
            }

            return result.second;
        }
                    
    

    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        bool LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        addDirichletBC (const dolfin::DirichletBC& dirichletCondition, 
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
            else
            {
                parameters ["system_is_assembled"] = false;
            }

            return result.second;
        }

    

    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        bool LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        addDirichletBC (dolfin::DirichletBC&& dirichletCondition, 
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
            else
            {
                parameters ["system_is_assembled"] = false;
            }

            return result.second;
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        bool LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        removeDirichletBC (const std::string& bcName)
        {
            dolfin::log (dolfin::DBG, "Removing dirichlet boundary condition \"%s\" from boundary conditions map...", 
                         bcName.c_str ());
            std::size_t nErasedElements = dirichletBCs_.erase (bcName);

            if (nErasedElements == 0)
            {
                dolfin::warning ("Cannot remove dirichlet boundary condition \"%s\" because it was not found in map",
                                 bcName.c_str ());
            }
            else
            {
                parameters ["system_is_assembled"] = false;
            }

            return nErasedElements == 1? true : false;
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        void LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        update ()
        {
            // define auxiliary string variables
            std::string desiredSolverType = parameters ["solver_type"];
            std::string desiredSolverMethod = parameters ["solver_method"];
            std::string desiredSolverPreconditioner = parameters ["solver_preconditioner"];
            
            bool needsSolverUpdating = (desiredSolverType != solverType_)
                                       || (desiredSolverMethod != solverMethod_) 
                                       || (desiredSolverPreconditioner != solverPreconditioner_);
            
            if (needsSolverUpdating)
            {
                dolfin::begin (dolfin::DBG, "Updating solver...");
                solver_ = createSolver_ ();
                dolfin::end ();
                
                parameters ["system_is_assembled"] = false;
            }
        }



    /******************* METHODS *******************/
    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        void LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        lumpMatrix () 
        {
            // create vector of 1s to get the sum of the rows by computing A*ones
            dolfin::Vector ones (MPI_COMM_WORLD, problemMatrix_.size (0));
            ones = 1;

            // vector to contain the result of the multiplication
            dolfin::Vector result;
            problemMatrix_.mult (ones, result);

            // set problemMatrix_ to diagonal
            problemMatrix_.zero ();
            problemMatrix_.set_diagonal (result);
        }
    


    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        void LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        assembleLinearSystem () 
        {
            // define auxiliary string variables
            bool systemIsAssembled = parameters ["system_is_assembled"];
            bool forceReassembleSystem = parameters ["force_reassemble_system"];
            bool needsReassembling = !systemIsAssembled || forceReassembleSystem;
            
            if (needsReassembling)
            {
                dolfin::log (dolfin::DBG, "Assembling bilinear form...");
                dolfin::assemble (problemMatrix_, bilinearForm_);
                
                dolfin::log (dolfin::DBG, "Assembling linear form...");
                dolfin::assemble (rhsVector_, linearForm_);
                
                if (!dirichletBCs_.empty ())
                {
                    dolfin::begin (dolfin::DBG, "Imposing Dirichlet's boundary conditions...");
                    for (auto &i : dirichletBCs_)
                    {
                        dolfin::log (dolfin::DBG, "Boundary condition: %s", i.first.c_str ());
                        i.second.apply (problemMatrix_, rhsVector_);
                    }
                    dolfin::end ();
                }

                if (bool(parameters["lump_matrix"]) == true)
                {
                    dolfin::begin (dolfin::DBG, "Lumping matrix...");
                    lumpMatrix ();
                    dolfin::end ();
                }
                
                parameters ["system_is_assembled"] = true;
            }
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        void LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        solve (const std::string& solveType) 
        {
            if (solveType != "default" && solveType != "stash")
            {
                dolfin::dolfin_error ("dcp: LinearProblem.h", 
                                      "solve",
                                      "Unknown solve type \"%s\" requested",
                                      solveType.c_str ());
            }
            
            dolfin::log (dolfin::DBG, "Solve type: %s", solveType.c_str ());
            
            update ();
            
            assembleLinearSystem ();
            
            dolfin::begin (dolfin::DBG, "Solving system...");
            solver_ -> set_operator (dolfin::reference_to_no_delete_pointer (problemMatrix_));
            if (solveType == "default")
            {
                solver_ -> solve (*(solution_.back ().second.vector ()), rhsVector_);
                
                // set stashedSolution_ to be equal to the last computed solution, so that when solution() is called
                // from a subiterations loop it gets the right one. In the case of "default" solveType, indeed, the
                // solution returned should be the same for all solution types
                stashedSolution_ = solution_.back ().second;
            }
            else if (solveType == "stash")
            {
                solver_ -> solve (*(stashedSolution_.vector ()), rhsVector_);
            }
            
            dolfin::end ();
        }
    
    
    
    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        dcp::LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>*
        LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        clone () const
        {
            dolfin::begin (dolfin::DBG, "Cloning object...");
            
            std::string cloneMethod = parameters["clone_method"];
            
            dolfin::log (dolfin::DBG, "Clone method: %s", cloneMethod.c_str ());
            dolfin::log (dolfin::DBG, "Creating new object of type LinearProblem...");
            
            // create new object
            dcp::LinearProblem <T_BilinearForm, T_LinearForm, T_LinearSolverFactory>* clonedProblem = nullptr;
            if (cloneMethod == "shallow_clone")
            {
                clonedProblem = 
                    new dcp::LinearProblem <T_BilinearForm, T_LinearForm, T_LinearSolverFactory> 
                    (this->functionSpace_,
                     this->bilinearForm_, 
                     this->linearForm_);
            }
            else if (cloneMethod == "deep_clone")
            {
                std::shared_ptr<dolfin::FunctionSpace> functionSpaceCopy (new dolfin::FunctionSpace (*functionSpace_));

                clonedProblem =
                    new dcp::LinearProblem <T_BilinearForm, T_LinearForm, T_LinearSolverFactory> 
                    (functionSpaceCopy,
                     this->bilinearForm_, 
                     this->linearForm_);
            }
            else
            {
                dolfin::dolfin_error ("dcp: LinearProblem.h",
                                      "clone",
                                      "Cannot clone linear differential problem. Unknown clone method: \"%s\"",
                                      cloneMethod.c_str ());
            }
            
            //copy dirichlet boundary conditions
            dolfin::log (dolfin::DBG, "Copying Dirichlet boundary conditions...");
            for (auto &i : this->dirichletBCs_)
            {
                clonedProblem->addDirichletBC (i.second, i.first);
            }

            // clear parameters set of newly created object so that it can be populated by the parameters of the object
            // being created. Set "system_is_assembled" to false, though, because in order for it to work the newly
            // created problem will have to reassemble matrix and rhs vector
            dolfin::log (dolfin::DBG, "Copying parameters to new object...");
            clonedProblem->parameters.clear ();
            clonedProblem->parameters = this->parameters;
            clonedProblem->update ();
            clonedProblem->parameters ["system_is_assembled"] = false;
            
            // copy solution
            dolfin::log (dolfin::DBG, "Copying solution...");
            clonedProblem->solution_ = this->solution_;
            
            dolfin::end ();
            
            return clonedProblem;
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        std::unique_ptr<dolfin::GenericLinearSolver> 
        LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        createSolver_ ()
        {
            std::string desiredSolverType = parameters["solver_type"];
            std::string desiredSolverMethod = parameters["solver_method"];
            std::string desiredSolverPreconditioner = parameters["solver_preconditioner"];
            
            std::unique_ptr<dolfin::GenericLinearSolver> solver;

            if (desiredSolverType == "lu_solver")
            {
                dolfin::log (dolfin::DBG, "Creating lu_solver...");
                solver.reset (new dolfin::LUSolver (desiredSolverMethod));
            }
            else if (desiredSolverType == "krylov_solver")
            {
                dolfin::log (dolfin::DBG, "Creating krylov_solver...");
                solver.reset (new dolfin::KrylovSolver (desiredSolverMethod, desiredSolverPreconditioner));
            }
            else
            {
                dolfin::log (dolfin::DBG, "Creating solver of type \"%s\"...", desiredSolverType.c_str ());
                dcp::LinearSolverFactory& factory = dcp::LinearSolverFactory::Instance ();
                solver = factory.create (desiredSolverType);
            }

            dolfin::log (dolfin::DBG, "Updating parameters...");

            // check if solver type has changed:
            // ----- if the solver is still the same, just update the parameters
            // ----- if the parameters do not have a parameters set corresponding to the one being created but
            //       current and desired solver set are the same, that means that there is no solver at all, so
            //       add solver parameters to parameters set
            // ----- otherwise, remove current solver parameters and add new ones
            if (parameters.has_parameter_set (desiredSolverType)) 
            {
                parameters.update (solver -> parameters);
                dolfin::log (dolfin::DBG, "Solver parameters have been updated.");
            }
            else if (solverType_ == desiredSolverType 
                     && solverMethod_ == desiredSolverMethod 
                     && solverPreconditioner_ == desiredSolverPreconditioner)
            {
                parameters.add (solver -> parameters);
                dolfin::log (dolfin::DBG, "Solver parameter set added to existing parameters.");
            }
            else
            {
                parameters.remove (solverType_);
                parameters.add (solver -> parameters);
                dolfin::log (dolfin::DBG, "Old solver parameter set replaced with new one.");
            }

            solverType_ = desiredSolverType;
            solverMethod_ = desiredSolverMethod;
            solverPreconditioner_ = desiredSolverPreconditioner;

                
            return solver;
        }
}
#endif
