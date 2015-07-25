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

#ifndef SRC_PROBLEMS_ALGEBRAICPROBLEM_H_INCLUDE_GUARD
#define SRC_PROBLEMS_ALGEBRAICPROBLEM_H_INCLUDE_GUARD

#include <dolfin/function/FunctionSpace.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/parameter/Parameters.h>
#include <vector>
#include <string>
#include <memory>
#include <dcp/problems/GenericProblem.h>
#include <dcp/expressions/GenericExpression.h>

namespace dcp
{
    /*! \class AlgebraicProblem AlgebraicProblem.h
     *  \brief Class for algebraic problems.
     *
     *  This class represents problem that are actually just operations on given coefficients.
     *  It inherits publicly from \c GenericProblem
     */

    class AlgebraicProblem : public dcp::GenericProblem
    {
        // ---------------------------------------------------------------------------------------------//  

        public:
            /******************* CONSTRUCTORS *******************/
            //! Default constructor is deleted. The class is not default constructable.
            AlgebraicProblem () = delete;

            //!  Constructor with shared pointers [1]
            /*!
             *  \param functionSpace the finite element space to be used to represent the solution
             *  \param expression the expression defining the problem. Its type must be a class derived from 
             *  \c dcp::GenericExpression
             *  
             *  The constructors also sets the following parameters:
             *      - \c "problem_type" a string describing the problem. Default value: \c "algebraic"
             */
            AlgebraicProblem (const std::shared_ptr<dolfin::FunctionSpace> functionSpace, 
                              const std::shared_ptr<dcp::GenericExpression> expression);


            //! Constructor with references [1]
            /*!
             *  \param functionSpace the finite element space to be used to represent the solution
             *  The stored function space's ownership will be unique to the object, since the protected member
             *  \param expression the expression defining the problem. Its type must be a class derived from 
             *  \c dcp::GenericExpression
             *  
             *  The constructors also sets the following parameters:
             *      - \c "problem_type" a string describing the problem. Default value: \c "algebraic"
             */
            AlgebraicProblem (const dolfin::FunctionSpace& functionSpace,
                              const dcp::GenericExpression& expression);


            /******************* DESTRUCTOR *******************/

            //! Destructor
            /*! 
             *  Default destructor, since members of the class are trivially 
             *  destructible.
             */
            virtual ~AlgebraicProblem () {};


            /******************* GETTERS *******************/
            //! Get const reference to the problem's expression
            /*! 
             *  \return a const reference to the problem's expression
             */
            virtual const dcp::GenericExpression& expression () const;

            /******************* SETTERS *******************/
            //! Set coefficient [1]. Override of virtual function in \c GenericProblem.
            /*!
             *  Only possible value for \c coefficientType is \c expression
             *  
             *  See \c GenericProblem documentation for more details on the function
             */
            virtual void setCoefficient (const std::string& coefficientType, 
                                         const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                         const std::string& coefficientName) override;

            //! Set coefficient [2]. Override of virtual function in \c GenericProblem.
            /*!
             *  Only possible value for \c coefficientType is \c expression
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


            /******************* METHODS *******************/
            //! Solve problem
            /*!
             *  This method solves the problem defined and stores the solution in the private member \c solution_. 
             *  See documentation of \c dcp::GenericProblem for more details on how the protected member \c solution_ 
             *  works and why it is declared as a \c std::pair.
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
            virtual dcp::AlgebraicProblem* clone () const override;

            // ---------------------------------------------------------------------------------------------//

        protected:
            //! The expression which defines the problem
            std::shared_ptr<dcp::GenericExpression> expression_;
            
            // ---------------------------------------------------------------------------------------------//

        private:
    };
}
#endif
