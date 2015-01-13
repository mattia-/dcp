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

#ifndef SRC_DIFFERENTIAL_PROBLEMS_NONLINEARPROBLEM_H_INCLUDE_GUARD
#define SRC_DIFFERENTIAL_PROBLEMS_NONLINEARPROBLEM_H_INCLUDE_GUARD

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/fem/solve.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/Vector.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/parameter/Parameters.h>
#include <vector>
#include <string>
#include <differential_problems/AbstractProblem.h>
#include <differential_problems/SubdomainType.h>

namespace dcp
{
    /*! \class NonlinearProblem NonlinearProblem.h
     *  \brief Class for non-linear differential problems.
     *
     *  This class represents problem of the form
     *  \f[
     *      \mbox{Find } u \in V : F \left(u, v\right) = 0 \ \forall\,v\,\in\,V
     *  \f]
     *  with \f$ F \left(u, v\right) : V \times V \rightarrow \mathds{R}\f$ non linear in the problem unknown \f$u\f$.
     *  
     *  It inherits publicly from \c AbstractProblem
     *  and it extends its functionalities to a concrete differential
     *  problem.
     *  Template arguments are:
     *  \arg T_ResidualForm the residual form type, that is the type that describes the form \f$F\f$
     *  \arg T_JacobianForm the jacobian form type, that is the type of the derivative of \f$F\f$ with respect to \f$u\f$,
     *  used in the Newton-Raphson method iterations
     */

    template <class T_ResidualForm_, class T_JacobianForm_>
        class NonlinearProblem : public dcp::AbstractProblem
        {
            // ---------------------------------------------------------------------------------------------//  

            public:
                /******************* TYPEDEFS *******************/
                typedef T_ResidualForm_ T_ResidualForm;
                typedef T_JacobianForm_ T_JacobianForm;
                
                
                /******************* CONSTRUCTORS *******************/
                //! Default constructor is deleted. The class is not default constructable.
                NonlinearProblem () = delete;

                //!  Constructor with shared pointers [1]
                /*!
                 *  \param mesh the problem mesh as a const \c std::shared_ptr to \c dolfin::Mesh
                 *  \param functionSpace the problem finite element space as a const \c std::shared_ptr to 
                 *  \c dolfin::FunctionSpace
                 *  \param residualFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the residual form 
                 *  \param jacobianFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the jacobian form. Default value is the empty string, in which case \c residualFormSolutionName
                 *  will be used for \c jacobianFormSolutionName
                 *  The stored mesh's and function space's ownership will be shared between the object and the input argument.
                 *  The residual and jacobian form will be created too, calling the constructor which takes the function space
                 *  as input.
                 */
                NonlinearProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
                                  const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
                                  const std::string& residualFormSolutionName,
                                  const std::string& jacobianFormSolutionName = "" );
                

                //! Constructor with references [1]
                /*!
                 *  \param mesh the problem mesh as a const \c dolfin::Mesh&
                 *  \param functionSpace the problem finite element space as a const \c dolfin::FunctionSpace&
                 *  \param residualFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the residual form 
                 *  \param jacobianFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the jacobian form. Default value is the empty string, in which case \c residualFormSolutionName
                 *  will be used for \c jacobianFormSolutionName
                 *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
                 *  initialized using the \c new operator and mesh's and functionSpace's copy constructor
                 *  The residual and jacobian form will be created too, calling the constructor which takes the function space
                 *  as input.
                 */
                NonlinearProblem (const dolfin::Mesh& mesh, 
                                  const dolfin::FunctionSpace& functionSpace,
                                  const std::string& residualFormSolutionName,
                                  const std::string& jacobianFormSolutionName = "" );

                //! Constructor with rvalue references [1]
                /*!
                 *  \param mesh the problem mesh as a dolfin::Mesh&&
                 *  \param functionSpace the problem finite element space as a dolfin::FunctionSpace&&
                 *  \param residualFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the residual form 
                 *  \param jacobianFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the jacobian form. Default value is the empty string, in which case \c residualFormSolutionName
                 *  will be used for \c jacobianFormSolutionName
                 *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
                 *  initialized using the \c new operator and mesh's and functionSpace's move constructor
                 *  The residual and jacobian form will be created too, calling the constructor which takes the function space
                 *  as input.
                 */
                NonlinearProblem (dolfin::Mesh&& mesh, 
                                  dolfin::FunctionSpace&& functionSpace,
                                  const std::string& residualFormSolutionName,
                                  const std::string& jacobianFormSolutionName = "" );

                
                //!  Constructor with shared pointers [2]
                /*!
                 *  \param mesh the problem mesh as a const \c std::shared_ptr to \c dolfin::Mesh
                 *  \param functionSpace the problem finite element space as a const \c std::shared_ptr 
                 *  \c to dolfin::FunctionSpace
                 *  \param residualForm a const reference to the problem's residual form
                 *  \param jacobianForm a const reference to the problem's jacobian form
                 *  \param residualFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the residual form 
                 *  \param jacobianFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the jacobian form. Default value is the empty string, in which case \c residualFormSolutionName
                 *  will be used for \c jacobianFormSolutionName
                 *  The stored mesh's and function space's ownership will be shared between the object and the input argument.
                 *  The residual and jacobian form will be created too, calling the constructor which takes the function space
                 *  as input.
                 */
                NonlinearProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
                                  const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
                                  const T_ResidualForm& residualForm,
                                  const T_JacobianForm& jacobianForm,
                                  const std::string& residualFormSolutionName,
                                  const std::string& jacobianFormSolutionName = "" );

                //! Constructor with references [2]
                /*!
                 *  \param mesh the problem mesh as a const \c dolfin::Mesh&
                 *  \param functionSpace the problem finite element space as a const \c dolfin::FunctionSpace&
                 *  \param residualForm a const reference to the problem's residual form
                 *  \param jacobianForm a const reference to the problem's jacobian form
                 *  \param residualFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the residual form 
                 *  \param jacobianFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the jacobian form. Default value is the empty string, in which case 
                 *  \c residualFormSolutionName will be used for \c jacobianFormSolutionName
                 *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
                 *  initialized using the \c new operator and mesh's and functionSpace's copy constructor
                 *  The residual and jacobian form will be created too, calling the constructor which takes the function space
                 *  as input.
                 */
                NonlinearProblem (const dolfin::Mesh& mesh, 
                                  const dolfin::FunctionSpace& functionSpace,
                                  const T_ResidualForm& residualForm,
                                  const T_JacobianForm& jacobianForm,
                                  const std::string& residualFormSolutionName,
                                  const std::string& jacobianFormSolutionName = "" );

                //! Constructor with rvalue references [2]
                /*!
                 *  \param mesh the problem mesh as a \c dolfin::Mesh&&
                 *  \param functionSpace the problem finite element space as a \c dolfin::FunctionSpace&&
                 *  \param residualForm a rvalue reference to the problem's residual form
                 *  \param jacobianForm a rvalue reference to the problem's jacobian form
                 *  \param residualFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the residual form
                 *  \param jacobianFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the jacobian form. Default value is the empty string, in which case 
                 *  \c residualFormSolutionName will be used for \c jacobianFormSolutionName
                 *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
                 *  initialized using the \c new operator and mesh's and functionSpace's move constructor
                 *  The residual and jacobian form will be created too, calling the constructor which takes the function space
                 *  as input.
                 */
                NonlinearProblem (dolfin::Mesh&& mesh, 
                                  dolfin::FunctionSpace&& functionSpace,
                                  T_ResidualForm&& residualForm,
                                  T_JacobianForm&& jacobianForm,
                                  const std::string& residualFormSolutionName,
                                  const std::string& jacobianFormSolutionName = "" );
                

                /******************* DESTRUCTOR *******************/
                
                //! Destructor
                /*! Default destructor, since members of the class are trivially 
                 * destructible.
                 */
                virtual ~NonlinearProblem () {};

                
                /******************* GETTERS *******************/
                //! Get const reference to the problem's residual form
                /*! 
                 *  \return a const reference to the problem's residual form
                 */
                virtual const T_ResidualForm& residualForm () const;

                //! Get const reference to the problem's jacobian form
                /*! 
                 *  \return a const reference to the problem's jacobian form
                 */
                virtual const T_JacobianForm& jacobianForm () const;

                /******************* SETTERS *******************/

                //! Set coefficient [1]. Override of virtual function in \c AbstractProblem.
                /*!
                 *  Possible values for \c coefficientType are:
                 *  \li \c residual_form to set the coefficient in the residual form
                 *  \li \c jacobian_form to set the coefficient in the jacobian form
                 *  \li \c initial_guess to set the initial guess for the nonlinear solver
                 *
                 *  Parameter \c coefficientName has default value set to \c default, which is only used when setting
                 *  initial guess and allows us to call the function as
                 *  \code
                 *  setCoefficient ("initial_guess", function);
                 *  \endcode
                 *  without specifying the coefficient name (which would not be used anyway).
                 *  Note that when setting the initial guess, even though parameter \c coefficientValue is of type
                 *  \c std::shared_ptr<const dolfin::GenericFunction>, the real type of the object pointed by 
                 *  \c coefficientValue can only be \c dolfin::Expression or \c dolfin::Function. This is because
                 *  \c setCoefficient() will call the assignement operator of class \c dolfin::Function, which only accepts
                 *  the two types mentioned before as input arguments.
                 *  See \c AbstractProblem documentation for more details on the function
                 */
                virtual void setCoefficient (const std::string& coefficientType, 
                                             const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                             const std::string& coefficientName = "default") override;

                //! Set coefficient [2]. Override of virtual function in \c AbstractProblem.
                /*!
                 *  Possible values for \c coefficientType are:
                 *  \li \c residual_form to set the coefficient in the residual form
                 *  \li \c jacobian_form to set the coefficient in the jacobian form
                 *  
                 *  See \c AbstractProblem documentation for more details on the function
                 */
                virtual void setCoefficient (const std::string& coefficientType,
                                             const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                             const std::size_t& coefficientNumber) override;

                //! Set integration subdomains for the forms. Override of virtual function in \c AbstractProblem
                /*! 
                 *  Possible values for \c formType are:
                 *  \li \c residual_form to set the integration subdomain in the residual form
                 *  \li \c jacobian_form to set the integration subdomain in the jacobian form
                 *  
                 *  See \c AbstractProblem documentation for more details on the function
                 */
                virtual void setIntegrationSubdomains (const std::string& formType,
                                                       std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                                       const dcp::SubdomainType& subdomainType) override;
                
                /******************* METHODS *******************/

                //! Solve problem
                /*!
                 * This method solves the problem defined. It uses the private members' value to set the problem and then
                 * stores the solution in the private member \c solution_
                 *
                 * \param type the solution type requested. In this class, the only possibility is to set
                 * \c type equal to \c "default".
                 */
                virtual void solve (const std::string& type = "default") override;

                //! Solve problem specifying flag
                /*!
                 *  \param solverParameters object of type \c dolfin::parameters that contain the parameters to
                 *  be used for the non linear solver
                 */
                virtual void solve (const dolfin::Parameters& solverParameters);

                //! Clone method. Overrides method in \c AbstractProblem
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
                virtual dcp::NonlinearProblem <T_ResidualForm, T_JacobianForm>* clone () const override;
                
                // ---------------------------------------------------------------------------------------------//

            protected:
                //! The residual form
                T_ResidualForm residualForm_;

                //! The jacobian form
                T_JacobianForm jacobianForm_;
                
                // ---------------------------------------------------------------------------------------------//

            private:
        };



    // ==============================================================================================//
    // ==================================== IMPLEMENTATION ==========================================//
    // ==============================================================================================//


    /******************* CONSTRUCTORS *******************/

    template <class T_ResidualForm, class T_JacobianForm>
        NonlinearProblem<T_ResidualForm, T_JacobianForm>::
        NonlinearProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
                          const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
                          const std::string& residualFormSolutionName,
                          const std::string& jacobianFormSolutionName) : 
            AbstractProblem (mesh, functionSpace),
            residualForm_ (*functionSpace_),
            jacobianForm_ (*functionSpace_, *functionSpace_)
        { 
            dolfin::begin (dolfin::DBG, "Building NonlinearProblem...");
            
            solution_.emplace_back (dolfin::Function (*functionSpace_));
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "nonlinear");
            parameters.add (dolfin::NonlinearVariationalSolver::default_parameters ());
            parameters.add ("residual_form_solution_name", residualFormSolutionName);
            if (jacobianFormSolutionName.empty ())
            {
                parameters.add ("jacobian_form_solution_name", residualFormSolutionName);
            }
            else
            {
                parameters.add ("jacobian_form_solution_name", jacobianFormSolutionName);
            }
            parameters.add ("solver_parameters_set_name", "nonlinear_variational_solver");
            parameters.add ("clone_method", "shallow_clone");
            
            dolfin::log (dolfin::DBG, "Setting initial guess...");
            if (jacobianFormSolutionName.empty ())
            {
                residualForm_.set_coefficient (residualFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
                jacobianForm_.set_coefficient (residualFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
            }
            else
            {
                residualForm_.set_coefficient (residualFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
                jacobianForm_.set_coefficient (jacobianFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
            }
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "NonlinearProblem object created");
        }



    template <class T_ResidualForm, class T_JacobianForm>
        NonlinearProblem<T_ResidualForm, T_JacobianForm>::
        NonlinearProblem (const dolfin::Mesh& mesh, 
                          const dolfin::FunctionSpace& functionSpace,
                          const std::string& residualFormSolutionName,
                          const std::string& jacobianFormSolutionName) : 
            AbstractProblem (mesh, functionSpace),
            residualForm_ (*functionSpace_),
            jacobianForm_ (*functionSpace_, *functionSpace_)
        { 
            dolfin::begin (dolfin::DBG, "Building NonlinearProblem...");
            
            solution_.emplace_back (dolfin::Function (*functionSpace_));
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "nonlinear");
            parameters.add (dolfin::NonlinearVariationalSolver::default_parameters ());
            parameters.add ("residual_form_solution_name", residualFormSolutionName);
            if (jacobianFormSolutionName.empty ())
            {
                parameters.add ("jacobian_form_solution_name", residualFormSolutionName);
            }
            else
            {
                parameters.add ("jacobian_form_solution_name", jacobianFormSolutionName);
            }
            parameters.add ("solver_parameters_set_name", "nonlinear_variational_solver");
            parameters.add ("clone_method", "shallow_clone");
            
            dolfin::log (dolfin::DBG, "Setting initial guess...");
            if (jacobianFormSolutionName.empty ())
            {
                residualForm_.set_coefficient (residualFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
                jacobianForm_.set_coefficient (residualFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
            }
            else
            {
                residualForm_.set_coefficient (residualFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
                jacobianForm_.set_coefficient (jacobianFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
            }
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "NonlinearProblem object created");
        }
            

    
    template <class T_ResidualForm, class T_JacobianForm>
        NonlinearProblem<T_ResidualForm, T_JacobianForm>::
        NonlinearProblem (dolfin::Mesh&& mesh, 
                          dolfin::FunctionSpace&& functionSpace,
                          const std::string& residualFormSolutionName,
                          const std::string& jacobianFormSolutionName) : 
            AbstractProblem (mesh, functionSpace),
            residualForm_ (*functionSpace_),
            jacobianForm_ (*functionSpace_, *functionSpace_)
        { 
            dolfin::begin (dolfin::DBG, "Building NonlinearProblem...");
            
            solution_.emplace_back (dolfin::Function (*functionSpace_));
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "nonlinear");
            parameters.add (dolfin::NonlinearVariationalSolver::default_parameters ());
            parameters.add ("residual_form_solution_name", residualFormSolutionName);
            if (jacobianFormSolutionName.empty ())
            {
                parameters.add ("jacobian_form_solution_name", residualFormSolutionName);
            }
            else
            {
                parameters.add ("jacobian_form_solution_name", jacobianFormSolutionName);
            }
            parameters.add ("solver_parameters_set_name", "nonlinear_variational_solver");
            parameters.add ("clone_method", "shallow_clone");
            
            dolfin::log (dolfin::DBG, "Setting initial guess...");
            if (jacobianFormSolutionName.empty ())
            {
                residualForm_.set_coefficient (residualFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
                jacobianForm_.set_coefficient (residualFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
            }
            else
            {
                residualForm_.set_coefficient (residualFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
                jacobianForm_.set_coefficient (jacobianFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
            }
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "NonlinearProblem object created");
        }
            


    template <class T_ResidualForm, class T_JacobianForm>
        NonlinearProblem<T_ResidualForm, T_JacobianForm>::
        NonlinearProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
                          const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
                          const T_ResidualForm& residualForm,
                          const T_JacobianForm& jacobianForm,
                          const std::string& residualFormSolutionName,
                          const std::string& jacobianFormSolutionName) : 
            AbstractProblem (mesh, functionSpace),
            residualForm_ (residualForm),
            jacobianForm_ (jacobianForm)
        { 
            dolfin::begin (dolfin::DBG, "Building NonlinearProblem...");
            
            solution_.emplace_back (dolfin::Function (*functionSpace_));
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "nonlinear");
            parameters.add (dolfin::NonlinearVariationalSolver::default_parameters ());
            parameters.add ("residual_form_solution_name", residualFormSolutionName);
            if (jacobianFormSolutionName.empty ())
            {
                parameters.add ("jacobian_form_solution_name", residualFormSolutionName);
            }
            else
            {
                parameters.add ("jacobian_form_solution_name", jacobianFormSolutionName);
            }
            parameters.add ("solver_parameters_set_name", "nonlinear_variational_solver");
            parameters.add ("clone_method", "shallow_clone");
            
            dolfin::log (dolfin::DBG, "Setting initial guess...");
            if (jacobianFormSolutionName.empty ())
            {
                residualForm_.set_coefficient (residualFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
                jacobianForm_.set_coefficient (residualFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
            }
            else
            {
                residualForm_.set_coefficient (residualFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
                jacobianForm_.set_coefficient (jacobianFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
            }
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "NonlinearProblem object created");
        }

    

    template <class T_ResidualForm, class T_JacobianForm>
        NonlinearProblem<T_ResidualForm, T_JacobianForm>::
        NonlinearProblem (const dolfin::Mesh& mesh, 
                          const dolfin::FunctionSpace& functionSpace,
                          const T_ResidualForm& residualForm,
                          const T_JacobianForm& jacobianForm,
                          const std::string& residualFormSolutionName,
                          const std::string& jacobianFormSolutionName) : 
            AbstractProblem (mesh, functionSpace),
            residualForm_ (residualForm),
            jacobianForm_ (jacobianForm)
        { 
            dolfin::begin (dolfin::DBG, "Building NonlinearProblem...");
            
            solution_.emplace_back (dolfin::Function (*functionSpace_));
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "nonlinear");
            parameters.add (dolfin::NonlinearVariationalSolver::default_parameters ());
            parameters.add ("residual_form_solution_name", residualFormSolutionName);
            if (jacobianFormSolutionName.empty ())
            {
                parameters.add ("jacobian_form_solution_name", residualFormSolutionName);
            }
            else
            {
                parameters.add ("jacobian_form_solution_name", jacobianFormSolutionName);
            }
            parameters.add ("solver_parameters_set_name", "nonlinear_variational_solver");
            parameters.add ("clone_method", "shallow_clone");
            
            dolfin::log (dolfin::DBG, "Setting initial guess...");
            if (jacobianFormSolutionName.empty ())
            {
                residualForm_.set_coefficient (residualFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
                jacobianForm_.set_coefficient (residualFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
            }
            else
            {
                residualForm_.set_coefficient (residualFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
                jacobianForm_.set_coefficient (jacobianFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
            }
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "NonlinearProblem object created");
        }
            


    template <class T_ResidualForm, class T_JacobianForm>
        NonlinearProblem<T_ResidualForm, T_JacobianForm>::
        NonlinearProblem (dolfin::Mesh&& mesh, 
                          dolfin::FunctionSpace&& functionSpace,
                          T_ResidualForm&& residualForm,
                          T_JacobianForm&& jacobianForm,
                          const std::string& residualFormSolutionName,
                          const std::string& jacobianFormSolutionName) : 
            AbstractProblem (mesh, functionSpace),
            residualForm_ (residualForm),
            jacobianForm_ (jacobianForm)
        { 
            dolfin::begin (dolfin::DBG, "Building NonlinearProblem...");
            
            solution_.emplace_back (dolfin::Function (*functionSpace_));
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "nonlinear");
            parameters.add (dolfin::NonlinearVariationalSolver::default_parameters ());
            parameters.add ("residual_form_solution_name", residualFormSolutionName);
            if (jacobianFormSolutionName.empty ())
            {
                parameters.add ("jacobian_form_solution_name", residualFormSolutionName);
            }
            else
            {
                parameters.add ("jacobian_form_solution_name", jacobianFormSolutionName);
            }
            parameters.add ("solver_parameters_set_name", "nonlinear_variational_solver");
            parameters.add ("clone_method", "shallow_clone");
            
            dolfin::log (dolfin::DBG, "Setting initial guess...");
            
            if (jacobianFormSolutionName.empty ())
            {
                residualForm_.set_coefficient (residualFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
                jacobianForm_.set_coefficient (residualFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
            }
            else
            {
                residualForm_.set_coefficient (residualFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
                jacobianForm_.set_coefficient (jacobianFormSolutionName, dolfin::reference_to_no_delete_pointer (solution_.back ()));
            }
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "NonlinearProblem object created");
        }
    


    /******************* GETTERS *******************/

    template <class T_ResidualForm, class T_JacobianForm>
        const T_ResidualForm& NonlinearProblem<T_ResidualForm, T_JacobianForm>::
        residualForm () const
        {
            return residualForm_;
        }



    template <class T_ResidualForm, class T_JacobianForm>
        const T_JacobianForm& NonlinearProblem<T_ResidualForm, T_JacobianForm>::
        jacobianForm () const
        {
            return jacobianForm_;
        }

    

    /******************* SETTERS *******************/

    template <class T_ResidualForm, class T_JacobianForm>
        void NonlinearProblem<T_ResidualForm, T_JacobianForm>::
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
            else if (coefficientType == "initial_guess")
            {
                dolfin::log (dolfin::DBG, "Setting initial guess...");
                // check whether coefficientValue is a pointer to dolfin::Function or dolfin::Expression
                if (std::dynamic_pointer_cast<const dolfin::Function> (coefficientValue) != nullptr)
                {
                    solution_.back () = *(std::dynamic_pointer_cast<const dolfin::Function> (coefficientValue));
                }
                else if (std::dynamic_pointer_cast<const dolfin::Expression> (coefficientValue) != nullptr)
                {
                    solution_.back () = *(std::dynamic_pointer_cast<const dolfin::Expression> (coefficientValue));
                }
                else
                {
                    dolfin::warning ("Cannot set initial guess in nonlinear differential problem. Input argument is \
                                     neither a dolfin::Function nor a dolfin::Expression");
                }
            }
            else
            {
                dolfin::warning ("Cannot set coefficient in non linear differential problem. Coefficient type \"%s\" unknown",
                                 coefficientType.c_str ());
            }
        }



    template <class T_ResidualForm, class T_JacobianForm>
        void NonlinearProblem<T_ResidualForm, T_JacobianForm>::
        setCoefficient (const std::string& coefficientType,
                        const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                        const std::size_t& coefficientNumber)
        {
            if (coefficientType == "residual_form")
            {
                dolfin::log (dolfin::DBG, "Setting residual form coefficient number %d...", coefficientNumber);
                residualForm_.set_coefficient (coefficientNumber, coefficientValue);
            }
            else if (coefficientType == "jacobian_form")
            {
                dolfin::log (dolfin::DBG, "Setting jacobian form coefficient number %d...", coefficientNumber);
                jacobianForm_.set_coefficient (coefficientNumber, coefficientValue);
            }
            else
            {
                dolfin::warning ("Cannot set coefficient in non linear differential problem. Coefficient type \"%s\" unknown",
                                 coefficientType.c_str ());
            }
        }



    template <class T_ResidualForm, class T_JacobianForm>
        void NonlinearProblem<T_ResidualForm, T_JacobianForm>::
        setIntegrationSubdomains (const std::string& formType,
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


    
    /******************* METHODS *******************/

    template <class T_ResidualForm, class T_JacobianForm>
        void NonlinearProblem<T_ResidualForm, T_JacobianForm>::
        solve (const std::string& type) 
        {
            if (type != "default")
            {
                dolfin::dolfin_error ("dcp: NonlinearProblem.h", 
                                      "solve",
                                      "Unknown solve type \"%s\" requested",
                                      type.c_str ());
            }
            
            dolfin::log (dolfin::DBG, "Solve type: %s", type.c_str ());
            
            dolfin::begin (dolfin::INFO, "Solving problem...");
            
            dolfin::log (dolfin::DBG, "Creating temporary vectors of dolfin::DirichletBC pointers...");
            // create vector of POINTERS to DirichletBC. This is needed to call the function dolfin::solve
            std::vector<const dolfin::DirichletBC*> tmpDirichletBCs (dirichletBCs_.size (), nullptr);
            std::size_t counter = 0;
            
            for (auto i = dirichletBCs_.begin (); i != dirichletBCs_.end (); ++i)
            {
                tmpDirichletBCs[counter] = &(i->second);
                counter++;
            }
            // note that when tmpDirichletBCs gets destroyd (upon exit of the function) the vector distructor will
            // call DirichletBC* destructor, which means that the pointers will be destroyed but the object they point
            // to will not
            
            // now solve non linear problem and store solution in solution_.
            // We need to check whether tmpDirichletBCs is empty, and call the appropriate function
            // Parameters passed to the solver are those stored in private member "parameters" in the set identified
            // by parameter "solver_parameters_set_name"
            std::string solverParametersSetName = parameters ["solver_parameters_set_name"];
            dolfin::log (dolfin::DBG, "Solver parameters set name is: %s", solverParametersSetName.c_str ());
            if (tmpDirichletBCs.size () != 0)
            {
                dolfin::solve (residualForm_ == 0, solution_.back (), tmpDirichletBCs, jacobianForm_, 
                               parameters (solverParametersSetName));
            }
            else
            {
                dolfin::solve (residualForm_ == 0, solution_.back (), jacobianForm_, 
                               parameters (solverParametersSetName));
            }
            
            dolfin::end ();
        }



    template <class T_ResidualForm, class T_JacobianForm>
        void NonlinearProblem<T_ResidualForm, T_JacobianForm>::
        solve (const dolfin::Parameters& solverParameters)
        {
            dolfin::begin (dolfin::DBG, "Solving problem...");
            
            dolfin::log (dolfin::DBG, "Creating temporary vectors of dolfin::DirichletBC pointers...");
            // create vector of POINTERS to DirichletBC. This is needed to call the function dolfin::solve
            std::vector<const dolfin::DirichletBC*> tmpDirichletBCs (dirichletBCs_.size (), nullptr);
            std::size_t counter = 0;
            
            for (auto i = dirichletBCs_.begin (); i != dirichletBCs_.end (); ++i)
            {
                tmpDirichletBCs[counter] = &(i->second);
                counter++;
            }
            // note that when tmpDirichletBCs gets destroyd (upon exit of the function) the vector distructor will
            // call DirichletBC* destructor, which means that the pointers will be destroyed but the object they point
            // to will not
            
            // now solve non linear problem and store solution in solution_.
            // We need to check whether tmpDirichletBCs is empty, and call the appropriate function
            // Parameters passed to the solver are those stored in the input variable
            if (tmpDirichletBCs.size () != 0)
            {
                dolfin::solve (residualForm_ == 0, solution_.back (), tmpDirichletBCs, jacobianForm_, solverParameters);
            }
            else
            {
                dolfin::solve (residualForm_ == 0, solution_.back (), jacobianForm_, solverParameters);
            }
            
            dolfin::end ();
        }
    
    
    
    template <class T_ResidualForm, class T_JacobianForm>
        dcp::NonlinearProblem <T_ResidualForm, T_JacobianForm>*
        NonlinearProblem<T_ResidualForm, T_JacobianForm>::
        clone () const
        {
            dolfin::begin (dolfin::DBG, "Cloning object...");
            
            std::string cloneMethod = parameters["clone_method"];
            
            dolfin::log (dolfin::DBG, "Clone method: %s", cloneMethod.c_str ());
            dolfin::log (dolfin::DBG, "Creating new object of type NonlinearProblem...");
            
            // create new object
            dcp::NonlinearProblem <T_ResidualForm, T_JacobianForm>* clonedProblem = nullptr;
            if (cloneMethod == "shallow_clone")
            {
                clonedProblem = 
                    new dcp::NonlinearProblem <T_ResidualForm, T_JacobianForm> 
                        ( this->mesh_,
                          this->functionSpace_,
                          this->residualForm_, 
                          this->jacobianForm_,
                         (this->parameters) ["residual_form_solution_name"],
                         (this->parameters) ["jacobian_form_solution_name"]
                        );
            }
            else if (cloneMethod == "deep_clone")
            {
                clonedProblem = 
                    new dcp::NonlinearProblem <T_ResidualForm, T_JacobianForm> 
                        (*(this->mesh_),
                         *(this->functionSpace_),
                           this->residualForm_, 
                           this->jacobianForm_,
                          (this->parameters) ["residual_form_solution_name"],
                          (this->parameters) ["jacobian_form_solution_name"]
                        );
            }
            else
            {
                dolfin::dolfin_error ("dcp: NonlinearProblem.h",
                                      "clone",
                                      "Cannot clone nonlinear differential problem. Unknown clone method: \"%s\"",
                                      cloneMethod.c_str ());
            }
            //copy dirichlet boundary conditions
            dolfin::log (dolfin::DBG, "Copying Dirichlet boundary conditions...");
            for (auto &i : this->dirichletBCs_)
            {
                clonedProblem->addDirichletBC (i.second, i.first);
            }
            
            // clear parameters set of newly created object so that it can be populated by the parameters of the object
            // being created
            dolfin::log (dolfin::DBG, "Copying parameters to new object...");
            clonedProblem->parameters.clear ();
            clonedProblem->parameters = this->parameters;
            
            // copy solution
            dolfin::log (dolfin::DBG, "Copying solution...");
            clonedProblem->solution_ = this->solution_;
            
            dolfin::end ();
            
            return clonedProblem;
        }
}
#endif

