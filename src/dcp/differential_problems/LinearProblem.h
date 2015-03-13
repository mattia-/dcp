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

#ifndef SRC_DIFFERENTIAL_PROBLEMS_LINEARPROBLEM_H_INCLUDE_GUARD
#define SRC_DIFFERENTIAL_PROBLEMS_LINEARPROBLEM_H_INCLUDE_GUARD

#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/la/LinearSolver.h>
#include <dolfin/la/GenericLinearSolver.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/Vector.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/parameter/Parameters.h>
#include <vector>
#include <string>
#include <memory>
#include <dcp/differential_problems/AbstractProblem.h>
#include <dcp/factories/LinearSolverFactory.h>
#include <dcp/differential_problems/SubdomainType.h>

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
     *  It inherits publicly from \c AbstractProblem
     *  and it extends its functionalities to a concrete differential
     *  problem.
     *  Template arguments are:
     *  \arg T_BilinearForm the bilinear form type
     *  \arg T_LinearForm the linear form type
     *  \arg T_LinearSolverFactory the type of the factory that creates the linear solver. By default, it is set
     *  to \c dcp::LinearSolverFactory
     */

    template <class T_BilinearForm_, class T_LinearForm_, class T_LinearSolverFactory_ = dcp::LinearSolverFactory>
        class LinearProblem : public dcp::AbstractProblem
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

                //!  Constructor with shared pointers [1]
                /*!
                 *  \param functionSpace the problem finite element space as a <tt> const std::shared_ptr </tt> to 
                 *  \c dolfin::FunctionSpace
                 *  The stored function space's ownership will be shared between the object and the input argument. The
                 *  bilinear and linear form will be created too, calling the constructor which takes the function space
                 *  as input.
                 *  The constructors also sets the following parameters:
                 *      - \c "problem_type" a string describing the problem. Default value: \c "linear"
                 *      - \c "current_solver_type" the type of the solver currently stored in the protected member 
                 *        \c solver_. Default value: \c "lu_solver"
                 *      - \c "current_solver_method" the method used by the solver currently stored in the protected 
                 *        member \c solver_. Default value: \c "default"
                 *      - \c "current_solver_preconditioner" the preconditioner used by the solver currently stored in 
                 *        the protected member \c solver_. Default value: \c "default"
                 *      - \c "desired_solver_type" the desired solver type. If it differs from 
                 *        \c "current_solver_type" when the \c solve method is called, the protected member \c solver_
                 *        will be updated accordingly. Default value: \c "lu_solver"
                 *      - \c "desired_solver_method" the desired solver method. If it differs from 
                 *        \c "current_solver_method" when the \c solve method is called, the protected member \c solver_
                 *        will be updated accordingly. Default value: \c "default"
                 *      - \c "desired_solver_preconditioner" the desired solver preconditioner. If it differs from 
                 *        \c "current_solver_preconditioner" when the \c solve method is called, the protected member 
                 *        \c solver_ will be updated accordingly. Default value: \c "default"
                 *      - \c "system_is_assembled" a flag that can be set to false if one wants to force a reassembly of
                 *        the linear system when the \ solve method is called. Default value: \c false
                 *      - \c "force_reassemble_system" a flag that, if set to \c true, causes the system to be 
                 *        reassembled every time the \c solve method is called. Default value: \c false
                 */
                LinearProblem (const std::shared_ptr<dolfin::FunctionSpace> functionSpace);


                //! Constructor with references [1]
                /*!
                 *  \param functionSpace the problem finite element space as a <tt> const dolfin::FunctionSpace& </tt>
                 *  The stored function space's ownership will be unique to the object, since the protected member
                 *  variable is initialized using the \c new operator and functionSpace's copy constructor. The
                 *  bilinear and linear form will be created too, calling the constructor which takes the function space
                 *  as input.
                 *  The constructors also sets the following parameters:
                 *      - \c "problem_type" a string describing the problem. Default value: \c "linear"
                 *      - \c "current_solver_type" the type of the solver currently stored in the protected member 
                 *        \c solver_. Default value: \c "lu_solver"
                 *      - \c "current_solver_method" the method used by the solver currently stored in the protected 
                 *        member \c solver_. Default value: \c "default"
                 *      - \c "current_solver_preconditioner" the preconditioner used by the solver currently stored in 
                 *        the protected member \c solver_. Default value: \c "default"
                 *      - \c "desired_solver_type" the desired solver type. If it differs from 
                 *        \c "current_solver_type" when the \c solve method is called, the protected member \c solver_
                 *        will be updated accordingly. Default value: \c "lu_solver"
                 *      - \c "desired_solver_method" the desired solver method. If it differs from 
                 *        \c "current_solver_method" when the \c solve method is called, the protected member \c solver_
                 *        will be updated accordingly. Default value: \c "default"
                 *      - \c "desired_solver_preconditioner" the desired solver preconditioner. If it differs from 
                 *        \c "current_solver_preconditioner" when the \c solve method is called, the protected member 
                 *        \c solver_ will be updated accordingly. Default value: \c "default"
                 *      - \c "system_is_assembled" a flag that can be set to false if one wants to force a reassembly of
                 *        the linear system when the \ solve method is called. Default value: \c false
                 *      - \c "force_reassemble_system" a flag that, if set to \c true, causes the system to be 
                 *        reassembled every time the \c solve method is called. Default value: \c false
                 */
                LinearProblem (const dolfin::FunctionSpace& functionSpace);

                //! Constructor with rvalue references [1]
                /*!
                 *  \param functionSpace the problem finite element space as a \c dolfin::FunctionSpace&&
                 *  The stored function space's ownership will be unique to the object, since the protected member
                 *  variable is initialized using the \c new operator and functionSpace's move constructor. The bilinear
                 *  and linear form will be created too, calling the constructor which takes the function space as
                 *  input.
                 *  The constructors also sets the following parameters:
                 *      - \c "problem_type" a string describing the problem. Default value: \c "linear"
                 *      - \c "current_solver_type" the type of the solver currently stored in the protected member 
                 *        \c solver_. Default value: \c "lu_solver"
                 *      - \c "current_solver_method" the method used by the solver currently stored in the protected 
                 *        member \c solver_. Default value: \c "default"
                 *      - \c "current_solver_preconditioner" the preconditioner used by the solver currently stored in 
                 *        the protected member \c solver_. Default value: \c "default"
                 *      - \c "desired_solver_type" the desired solver type. If it differs from 
                 *        \c "current_solver_type" when the \c solve method is called, the protected member \c solver_
                 *        will be updated accordingly. Default value: \c "lu_solver"
                 *      - \c "desired_solver_method" the desired solver method. If it differs from 
                 *        \c "current_solver_method" when the \c solve method is called, the protected member \c solver_
                 *        will be updated accordingly. Default value: \c "default"
                 *      - \c "desired_solver_preconditioner" the desired solver preconditioner. If it differs from 
                 *        \c "current_solver_preconditioner" when the \c solve method is called, the protected member 
                 *        \c solver_ will be updated accordingly. Default value: \c "default"
                 *      - \c "system_is_assembled" a flag that can be set to false if one wants to force a reassembly of
                 *        the linear system when the \ solve method is called. Default value: \c false
                 *      - \c "force_reassemble_system" a flag that, if set to \c true, causes the system to be 
                 *        reassembled every time the \c solve method is called. Default value: \c false
                 */
                LinearProblem (dolfin::FunctionSpace&& functionSpace);

                //!  Constructor with shared pointers [2]
                /*!
                 *  \param functionSpace the problem finite element space as a <tt>const std::shared_ptr </tt> to 
                 *  \c dolfin::FunctionSpace
                 *  \param bilinearForm a \c const reference to the problem's bilinear form
                 *  \param linearForm a \c const reference to the problem's linear form
                 *  The stored function space's ownership will be shared between the object and the input argument.  The
                 *  bilinear and linear form will be created too, calling the constructor which takes the function space
                 *  as input.
                 *  The constructors also sets the following parameters:
                 *      - \c "problem_type" a string describing the problem. Default value: \c "linear"
                 *      - \c "current_solver_type" the type of the solver currently stored in the protected member 
                 *        \c solver_. Default value: \c "lu_solver"
                 *      - \c "current_solver_method" the method used by the solver currently stored in the protected 
                 *        member \c solver_. Default value: \c "default"
                 *      - \c "current_solver_preconditioner" the preconditioner used by the solver currently stored in 
                 *        the protected member \c solver_. Default value: \c "default"
                 *      - \c "desired_solver_type" the desired solver type. If it differs from 
                 *        \c "current_solver_type" when the \c solve method is called, the protected member \c solver_
                 *        will be updated accordingly. Default value: \c "lu_solver"
                 *      - \c "desired_solver_method" the desired solver method. If it differs from 
                 *        \c "current_solver_method" when the \c solve method is called, the protected member \c solver_
                 *        will be updated accordingly. Default value: \c "default"
                 *      - \c "desired_solver_preconditioner" the desired solver preconditioner. If it differs from 
                 *        \c "current_solver_preconditioner" when the \c solve method is called, the protected member 
                 *        \c solver_ will be updated accordingly. Default value: \c "default"
                 *      - \c "system_is_assembled" a flag that can be set to false if one wants to force a reassembly of
                 *        the linear system when the \ solve method is called. Default value: \c false
                 *      - \c "force_reassemble_system" a flag that, if set to \c true, causes the system to be 
                 *        reassembled every time the \c solve method is called. Default value: \c false
                 */
                LinearProblem (const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
                               const T_BilinearForm& bilinearForm,
                               const T_LinearForm& linearForm);

                //! Constructor with references [2]
                /*!
                 *  \param functionSpace the problem finite element space as a <tt>const dolfin::FunctionSpace& </tt>
                 *  \param bilinearForm a \c const reference to the problem's bilinear form
                 *  \param linearForm a \c const reference to the problem's linear form
                 *  The stored function space's ownership will be unique to the object, since the protected member
                 *  variable is initialized using the \c new operator and functionSpace's copy constructor. The bilinear
                 *  and linear form will be created too, calling the constructor which takes the function space as
                 *  input.
                 *  The constructors also sets the following parameters:
                 *      - \c "problem_type" a string describing the problem. Default value: \c "linear"
                 *      - \c "current_solver_type" the type of the solver currently stored in the protected member 
                 *        \c solver_. Default value: \c "lu_solver"
                 *      - \c "current_solver_method" the method used by the solver currently stored in the protected 
                 *        member \c solver_. Default value: \c "default"
                 *      - \c "current_solver_preconditioner" the preconditioner used by the solver currently stored in 
                 *        the protected member \c solver_. Default value: \c "default"
                 *      - \c "desired_solver_type" the desired solver type. If it differs from 
                 *        \c "current_solver_type" when the \c solve method is called, the protected member \c solver_
                 *        will be updated accordingly. Default value: \c "lu_solver"
                 *      - \c "desired_solver_method" the desired solver method. If it differs from 
                 *        \c "current_solver_method" when the \c solve method is called, the protected member \c solver_
                 *        will be updated accordingly. Default value: \c "default"
                 *      - \c "desired_solver_preconditioner" the desired solver preconditioner. If it differs from 
                 *        \c "current_solver_preconditioner" when the \c solve method is called, the protected member 
                 *        \c solver_ will be updated accordingly. Default value: \c "default"
                 *      - \c "system_is_assembled" a flag that can be set to false if one wants to force a reassembly of
                 *        the linear system when the \ solve method is called. Default value: \c false
                 *      - \c "force_reassemble_system" a flag that, if set to \c true, causes the system to be 
                 *        reassembled every time the \c solve method is called. Default value: \c false
                 */
                LinearProblem (const dolfin::FunctionSpace& functionSpace,
                               const T_BilinearForm& bilinearForm,
                               const T_LinearForm& linearForm);

                //! Constructor with rvalue references [2]
                /*!
                 *  \param functionSpace the problem finite element space as a \c dolfin::FunctionSpace&&
                 *  \param bilinearForm a rvalue reference to the problem's bilinear form
                 *  \param linearForm a rvalue reference to the problem's linear form
                 *  The stored function space's ownership will be unique to the object, since the protected member
                 *  variable is initialized using the \c new operator and functionSpace's move constructor The bilinear
                 *  and linear form will be created too, calling the constructor which takes the function space as
                 *  input.
                 *  The constructors also sets the following parameters:
                 *      - \c "problem_type" a string describing the problem. Default value: \c "linear"
                 *      - \c "current_solver_type" the type of the solver currently stored in the protected member 
                 *        \c solver_. Default value: \c "lu_solver"
                 *      - \c "current_solver_method" the method used by the solver currently stored in the protected 
                 *        member \c solver_. Default value: \c "default"
                 *      - \c "current_solver_preconditioner" the preconditioner used by the solver currently stored in 
                 *        the protected member \c solver_. Default value: \c "default"
                 *      - \c "desired_solver_type" the desired solver type. If it differs from 
                 *        \c "current_solver_type" when the \c solve method is called, the protected member \c solver_
                 *        will be updated accordingly. Default value: \c "lu_solver"
                 *      - \c "desired_solver_method" the desired solver method. If it differs from 
                 *        \c "current_solver_method" when the \c solve method is called, the protected member \c solver_
                 *        will be updated accordingly. Default value: \c "default"
                 *      - \c "desired_solver_preconditioner" the desired solver preconditioner. If it differs from 
                 *        \c "current_solver_preconditioner" when the \c solve method is called, the protected member 
                 *        \c solver_ will be updated accordingly. Default value: \c "default"
                 *      - \c "system_is_assembled" a flag that can be set to false if one wants to force a reassembly of
                 *        the linear system when the \ solve method is called. Default value: \c false
                 *      - \c "force_reassemble_system" a flag that, if set to \c true, causes the system to be 
                 *        reassembled every time the \c solve method is called. Default value: \c false
                 */
                LinearProblem (dolfin::FunctionSpace&& functionSpace,
                               T_BilinearForm&& bilinearForm,
                               T_LinearForm&& linearForm);

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
                virtual const T_BilinearForm& bilinearForm () const;

                //! Get const reference to the problem's linear form
                /*! 
                 *  \return a const reference to the problem's linear form
                 */
                virtual const T_LinearForm& linearForm () const;

                //! Get const reference to the problem's linear operator
                /*!
                 *  \return a const reference to the problem's linear operator, which
                 *  is a \c dolfin::Matrix
                 */
                virtual const dolfin::Matrix& linearOperator () const;

                //! Get const reference to the problem's right hand side
                /*!
                 *  \return a const reference to the problem's right hand side, which
                 *  is a \c dolfin::Vector
                 */
                virtual const dolfin::Vector& rhs () const;


                /******************* SETTERS *******************/
                
                //! Set coefficient [1]. Override of virtual function in \c AbstractProblem.
                /*!
                 *  Possible values for \c coefficientType are:
                 *  \li \c bilinear_form to set the coefficient in the bilinear form
                 *  \li \c linear_form to set the coefficient in the linear form
                 *  
                 *  See \c AbstractProblem documentation for more details on the function
                 */
                virtual void setCoefficient (const std::string& coefficientType, 
                                             const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                             const std::string& coefficientName) override;

                //! Set coefficient [2]. Override of virtual function in \c AbstractProblem.
                /*!
                 *  Possible values for \c coefficientType are:
                 *  \li \c bilinear_form to set the coefficient in the bilinear form
                 *  \li \c linear_form to set the coefficient in the linear form
                 *  
                 *  See \c AbstractProblem documentation for more details on the function
                 */
                virtual void setCoefficient (const std::string& coefficientType,
                                             const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                             const std::size_t& coefficientNumber) override;

                //! Set integration subdomains for the forms. Override of virtual function in \c AbstractProblem
                /*! 
                 *  Possible values for \c formType are:
                 *  \li \c bilinear_form to set the integration subdomain in the bilinear form
                 *  \li \c linear_form to set the integration subdomain in the linear form
                 *  
                 *  See \c AbstractProblem documentation for more details on the function
                 */
                virtual void setIntegrationSubdomain (const std::string& formType,
                                                       std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                                       const dcp::SubdomainType& subdomainType) override;

                //! Add Dirichlet boundary condition to the problem [1]. Overrides method in \c AbstractProblem
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

                //! Add Dirichlet boundary condition to the problem [2]. Overrides method in \c AbstractProblem
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

                //! Add Dirichlet boundary condition to the problem [5]. Overrides method in \c AbstractProblem
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

                //! Add Dirichlet boundary condition to the problem [6]. Overrides method in \c AbstractProblem
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

                //! Remove Dirichlet boundary condition with given position. Overrides method in \c AbstractProblem
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

                //! Solve problem
                /*!
                 *  This method solves the problem defined. It uses the private members' value to set the problem and then
                 *  stores the solution in the private member \c solution_. See documentation of \c dcp::AbstractProblem
                 *  for more details on how the protected member \c solution_ works and why it is declared as a 
                 *  \c std::pair.
                 *  Note that the method \c solve() checks the problem's member
                 *  \c parameters to decided whether problem's matrix and vector should be reassembled and if
                 *  the values of the parameters \c desired_solver_type, \c desired_solver_method and 
                 *  \c desired_solver_preconditioner match the values of \c current_solver_type, \c current_solver_method
                 *  and \c current_solver_preconditioner. If they differ, it calls \c createSolver().
                 * 
                 *  \param type the solution type wanted. Possible values are:
                 *  \li \c "default" : the normal solution process
                 *  \li \c "must_reassemble" : force system reassembly
                 */
                virtual void solve (const std::string& type = "default") override;

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
                virtual dcp::LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>* clone () const override;

                // ---------------------------------------------------------------------------------------------//

            protected:
                //! Creates a linear solver object of type passed as input. 
                /*!
                 *  The solver will be created using the parameters \c desired_solver_type, \c desired_solver_method
                 *  and \c desired_solver_preconditioner set in the protected member \c parameters. It will
                 *  also set the parameters \c current_solver_type, \c current_solver_method and 
                 *  \c current_solver_preconditioner in the same set of parameters substituting the current values with 
                 *  the values being used to create the solver
                 */
                std::unique_ptr<dolfin::GenericLinearSolver> createSolver ();
                
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
                /*!
                 *  We use a \c std::shared_ptr for compatibility with FEniCS, version 1.3.0.
                 *  Note that there is no write access to this variable from outside the class,
                 *  so it is guaranteed to be a unique shared pointer
                 */
                std::shared_ptr<dolfin::Matrix> problemMatrix_;

                //! Vector to store the right hand side of the discrete problem
                dolfin::Vector rhsVector_;

                // ---------------------------------------------------------------------------------------------//

            private:
        };



    // ============================================================================================== //
    // ==================================== IMPLEMENTATION ========================================== //
    // ============================================================================================== //

    /******************* CONSTRUCTORS *******************/

    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        LinearProblem (const std::shared_ptr<dolfin::FunctionSpace> functionSpace) :
            AbstractProblem (functionSpace),
            bilinearForm_ (*functionSpace_, *functionSpace_),
            linearForm_ (*functionSpace_),
            solver_ (nullptr),
            problemMatrix_ (new dolfin::Matrix),
            rhsVector_ ()
        { 
            std::string solverType = "lu_solver";
            std::string solverMethod = "default";
            std::string solverPreconditioner = "default";

            dolfin::begin (dolfin::DBG, "Building LinearProblem...");

            solution_.emplace_back (std::make_pair (-1, dolfin::Function (*functionSpace_)));

            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "linear");
            parameters.add ("current_solver_type", solverType);
            parameters.add ("current_solver_method", solverMethod);
            parameters.add ("current_solver_preconditioner", solverPreconditioner);
            parameters.add ("desired_solver_type", solverType);
            parameters.add ("desired_solver_method", solverMethod);
            parameters.add ("desired_solver_preconditioner", solverPreconditioner);
            parameters.add ("system_is_assembled", false);
            parameters.add ("force_reassemble_system", false);  

            dolfin::begin (dolfin::DBG, "Creating solver...");
            solver_ = createSolver ();
            dolfin::end ();

            dolfin::end ();

            dolfin::log (dolfin::DBG, "LinearProblem object created");
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        LinearProblem (const dolfin::FunctionSpace& functionSpace) :
            AbstractProblem (functionSpace),
            bilinearForm_ (*functionSpace_, *functionSpace_),
            linearForm_ (*functionSpace_),
            solver_ (nullptr),
            problemMatrix_ (new dolfin::Matrix),
            rhsVector_ ()
        { 
            std::string solverType = "lu_solver";
            std::string solverMethod = "default";
            std::string solverPreconditioner = "default";

            dolfin::begin (dolfin::DBG, "Building LinearProblem...");
            
            solution_.emplace_back (std::make_pair (-1, dolfin::Function (*functionSpace_)));
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "linear");
            parameters.add ("current_solver_type", solverType);
            parameters.add ("current_solver_method", solverMethod);
            parameters.add ("current_solver_preconditioner", solverPreconditioner);
            parameters.add ("desired_solver_type", solverType);
            parameters.add ("desired_solver_method", solverMethod);
            parameters.add ("desired_solver_preconditioner", solverPreconditioner);
            parameters.add ("system_is_assembled", false);
            parameters.add ("force_reassemble_system", false);  
            
            dolfin::begin (dolfin::DBG, "Creating solver...");
            solver_ = createSolver ();
            dolfin::end ();
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "LinearProblem object created");
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        LinearProblem (dolfin::FunctionSpace&& functionSpace) :
            AbstractProblem (functionSpace),
            bilinearForm_ (*functionSpace_, *functionSpace_),
            linearForm_ (*functionSpace_),
            solver_ (nullptr),
            problemMatrix_ (new dolfin::Matrix),
            rhsVector_ ()
        { 
            std::string solverType = "lu_solver";
            std::string solverMethod = "default";
            std::string solverPreconditioner = "default";

            dolfin::begin (dolfin::DBG, "Building LinearProblem...");
            
            solution_.emplace_back (std::make_pair (-1, dolfin::Function (*functionSpace_)));
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "linear");
            parameters.add ("current_solver_type", solverType);
            parameters.add ("current_solver_method", solverMethod);
            parameters.add ("current_solver_preconditioner", solverPreconditioner);
            parameters.add ("desired_solver_type", solverType);
            parameters.add ("desired_solver_method", solverMethod);
            parameters.add ("desired_solver_preconditioner", solverPreconditioner);
            parameters.add ("system_is_assembled", false);
            parameters.add ("force_reassemble_system", false);
        
            dolfin::begin (dolfin::DBG, "Creating solver...");
            solver_ = createSolver ();
            dolfin::end ();
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "LinearProblem object created");
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        LinearProblem (const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
                       const T_BilinearForm& bilinearForm,
                       const T_LinearForm& linearForm) :
            AbstractProblem (functionSpace),
            bilinearForm_ (bilinearForm),
            linearForm_ (linearForm),
            solver_ (nullptr),
            problemMatrix_ (new dolfin::Matrix),
            rhsVector_ ()
        { 
            std::string solverType = "lu_solver";
            std::string solverMethod = "default";
            std::string solverPreconditioner = "default";

            dolfin::begin (dolfin::DBG, "Building LinearProblem...");
            
            solution_.emplace_back (std::make_pair (-1, dolfin::Function (*functionSpace_)));
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "linear");
            parameters.add ("current_solver_type", solverType);
            parameters.add ("current_solver_method", solverMethod);
            parameters.add ("current_solver_preconditioner", solverPreconditioner);
            parameters.add ("desired_solver_type", solverType);
            parameters.add ("desired_solver_method", solverMethod);
            parameters.add ("desired_solver_preconditioner", solverPreconditioner);
            parameters.add ("system_is_assembled", false);
            parameters.add ("force_reassemble_system", false);
            
            dolfin::begin (dolfin::DBG, "Creating solver...");
            solver_ = createSolver ();
            dolfin::end ();
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "LinearProblem object created");
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        LinearProblem (const dolfin::FunctionSpace& functionSpace,
                       const T_BilinearForm& bilinearForm,
                       const T_LinearForm& linearForm) :
            AbstractProblem (functionSpace),
            bilinearForm_ (bilinearForm),
            linearForm_ (linearForm),
            solver_ (nullptr),
            problemMatrix_ (new dolfin::Matrix),
            rhsVector_ ()
        { 
            std::string solverType = "lu_solver";
            std::string solverMethod = "default";
            std::string solverPreconditioner = "default";

            dolfin::begin (dolfin::DBG, "Building LinearProblem...");
            
            solution_.emplace_back (std::make_pair (-1, dolfin::Function (*functionSpace_)));
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "linear");
            parameters.add ("current_solver_type", solverType);
            parameters.add ("current_solver_method", solverMethod);
            parameters.add ("current_solver_preconditioner", solverPreconditioner);
            parameters.add ("desired_solver_type", solverType);
            parameters.add ("desired_solver_method", solverMethod);
            parameters.add ("desired_solver_preconditioner", solverPreconditioner);
            parameters.add ("system_is_assembled", false);
            parameters.add ("force_reassemble_system", false);  
            
            dolfin::begin (dolfin::DBG, "Creating solver...");
            solver_ = createSolver ();
            dolfin::end ();
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "LinearProblem object created");
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        LinearProblem (dolfin::FunctionSpace&& functionSpace,
                       T_BilinearForm&& bilinearForm,
                       T_LinearForm&& linearForm) :
            AbstractProblem (functionSpace),
            bilinearForm_ (std::move (bilinearForm)),
            linearForm_ (std::move (linearForm)),
            solver_ (nullptr),
            problemMatrix_ (new dolfin::Matrix),
            rhsVector_ ()
        {
            std::string solverType = "lu_solver";
            std::string solverMethod = "default";
            std::string solverPreconditioner = "default";

            dolfin::begin (dolfin::DBG, "Building LinearProblem...");
            
            solution_.emplace_back (std::make_pair (-1, dolfin::Function (*functionSpace_)));
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "linear");
            parameters.add ("current_solver_type", solverType);
            parameters.add ("current_solver_method", solverMethod);
            parameters.add ("current_solver_preconditioner", solverPreconditioner);
            parameters.add ("desired_solver_type", solverType);
            parameters.add ("desired_solver_method", solverMethod);
            parameters.add ("desired_solver_preconditioner", solverPreconditioner);
            parameters.add ("system_is_assembled", false);
            parameters.add ("force_reassemble_system", false);  
            
            dolfin::begin (dolfin::DBG, "Creating solver...");
            solver_ = createSolver ();
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
            return *problemMatrix_;
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
                dolfin::warning ("Cannot set coefficient in linear differential problem. Form type \"%s\" unknown",
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
                dolfin::warning ("Cannot set coefficient in linear differential problem. Form type \"%s\" unknown",
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
                dolfin::warning ("Dirichlet boundary condition \"%s\" not found in map", bcName.c_str ());
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
            std::string desiredSolverType = parameters ["desired_solver_type"];
            std::string desiredSolverMethod = parameters ["desired_solver_method"];
            std::string desiredSolverPreconditioner = parameters ["desired_solver_preconditioner"];
            
            std::string currentSolverType = parameters ["current_solver_type"];
            std::string currentSolverMethod = parameters ["current_solver_method"];
            std::string currentSolverPreconditioner = parameters ["current_solver_preconditioner"];
            
            bool needsSolverUpdating = (desiredSolverType != currentSolverType)
                                       || (desiredSolverMethod != currentSolverMethod) 
                                       || (desiredSolverPreconditioner != currentSolverPreconditioner);
            
            if (needsSolverUpdating)
            {
                dolfin::begin (dolfin::DBG, "Updating solver...");
                solver_ = createSolver ();
                dolfin::end ();
                
                parameters ["system_is_assembled"] = false;
            }
        }



    /******************* METHODS *******************/
    
    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        void LinearProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        solve (const std::string& type) 
        {
            if (type != "default" && type != "must_reassemble")
            {
                dolfin::dolfin_error ("dcp: LinearProblem.h", 
                                      "solve",
                                      "Unknown solve type \"%s\" requested",
                                      type.c_str ());
            }
            
            dolfin::log (dolfin::DBG, "Solve type: %s", type.c_str ());
            
            update ();
            
            dolfin::begin (dolfin::INFO, "Solving problem...");
            
            // variable to save current value of parameter force_reassemble_system
            bool forceReassembleSystemBackup = parameters ["force_reassemble_system"];
            if (type == "must_reassemble")
            {
                // overwrite parameter force_reassemble_system with input value and solve problem
                parameters ["force_reassemble_system"] = true;
            }
            
            // define auxiliary string variables
            bool systemIsAssembled = parameters ["system_is_assembled"];
            bool forceReassembleSystem = parameters ["force_reassemble_system"];
            bool needsReassembling = !systemIsAssembled || forceReassembleSystem;
            
            if (needsReassembling)
            {
                dolfin::begin (dolfin::DBG, "Assembling system...");
                
                dolfin::log (dolfin::DBG, "Assembling bilinear form...");
                dolfin::assemble (*problemMatrix_, bilinearForm_);
                
                dolfin::log (dolfin::DBG, "Assembling linear form...");
                dolfin::assemble (rhsVector_, linearForm_);
                
                if (!dirichletBCs_.empty ())
                {
                    dolfin::begin (dolfin::DBG, "Imposing Dirichlet's boundary conditions...");
                    for (auto &i : dirichletBCs_)
                    {
                        dolfin::begin (dolfin::DBG, "Boundary condition: %s", i.first.c_str ());
                        i.second.apply (*problemMatrix_, rhsVector_);
                        dolfin::end ();
                    }
                    dolfin::end ();
                }
                
                solver_ -> set_operator (problemMatrix_);
                
                parameters ["system_is_assembled"] = true;
                
                dolfin::end ();
            }
            
            solver_ -> solve (*(solution_.back ().second.vector ()), rhsVector_);
            
            if (type == "must_reassemble")
            {
                // restore value of parameter force_reassemble_system
                parameters ["force_reassemble_system"] = forceReassembleSystemBackup;
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
                clonedProblem =
                    new dcp::LinearProblem <T_BilinearForm, T_LinearForm, T_LinearSolverFactory> 
                    (*(this->functionSpace_),
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
        createSolver ()
        {
            std::string currentSolverType = parameters["current_solver_type"];
            std::string currentSolverMethod = parameters["current_solver_method"];
            std::string currentSolverPreconditioner = parameters["current_solver_preconditioner"];
            std::string desiredSolverType = parameters["desired_solver_type"];
            std::string desiredSolverMethod = parameters["desired_solver_method"];
            std::string desiredSolverPreconditioner = parameters["desired_solver_preconditioner"];
            
            if (desiredSolverType == "lu_solver")
            {
                dolfin::log (dolfin::DBG, "Creating lu_solver...");
                std::unique_ptr<dolfin::GenericLinearSolver> solver (new dolfin::LUSolver (desiredSolverMethod));
                
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
                }
                else if (currentSolverType == desiredSolverType 
                         && currentSolverMethod == desiredSolverMethod 
                         && currentSolverPreconditioner == desiredSolverPreconditioner)
                {
                    parameters.add (solver -> parameters);
                }
                else
                {
                    parameters.remove (currentSolverType);
                    parameters.add (solver -> parameters);
                }
                
                parameters ["current_solver_type"] = desiredSolverType;
                parameters ["current_solver_method"] = desiredSolverMethod;
                parameters ["current_solver_preconditioner"] = desiredSolverPreconditioner;

                
                return solver;
            }
            else if (desiredSolverType == "krylov_solver")
            {
                dolfin::log (dolfin::DBG, "Creating krylov_solver...");
                std::unique_ptr<dolfin::GenericLinearSolver> solver (new dolfin::KrylovSolver (desiredSolverMethod, 
                                                                                               desiredSolverPreconditioner));
                
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
                }
                else if (currentSolverType == desiredSolverType 
                         && currentSolverMethod == desiredSolverMethod 
                         && currentSolverPreconditioner == desiredSolverPreconditioner)
                {
                    parameters.add (solver -> parameters);
                }
                else
                {
                    parameters.remove (currentSolverType);
                    parameters.add (solver -> parameters);
                }
                
                parameters ["current_solver_type"] = desiredSolverType;
                parameters ["current_solver_method"] = desiredSolverMethod;
                parameters ["current_solver_preconditioner"] = desiredSolverPreconditioner;

                
                return solver;
            }
            else
            {
                dolfin::log (dolfin::DBG, "Creating solver of type \"%s\"...", desiredSolverType.c_str ());
                dcp::LinearSolverFactory& factory = dcp::LinearSolverFactory::Instance ();
                auto solver = factory.create (desiredSolverType);
                
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
                }
                else if (currentSolverType == desiredSolverType 
                         && currentSolverMethod == desiredSolverMethod 
                         && currentSolverPreconditioner == desiredSolverPreconditioner)
                {
                    parameters.add (solver -> parameters);
                }
                else
                {
                    parameters.remove (currentSolverType);
                    parameters.add (solver -> parameters);
                }
                
                parameters ["current_solver_type"] = desiredSolverType;
                parameters ["current_solver_method"] = desiredSolverMethod;
                parameters ["current_solver_preconditioner"] = desiredSolverPreconditioner;

                
                return solver;
            }
        }
}
#endif
