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

#ifndef SRC_SPLITTING_METHODS_ABSTRACTSPLITTINGMETHOD_H_INCLUDE_GUARD
#define SRC_SPLITTING_METHODS_ABSTRACTSPLITTINGMETHOD_H_INCLUDE_GUARD

#include <vector>
#include <memory>
#include <initializer_list>
#include <string>
#include <dcp/differential_problems/SubdomainType.h>
#include <dcp/differential_problems/TimeDependentEquationSystem.h>


namespace dcp
{
    /*! \class AbstractSplittingMethod AbstractSplittingMethod.h
     *  \brief Abstract base class for splitting methods for a system of 
     *  differential problems. 
     *         
     *  This class contains the basic interface for a splitting method to 
     *  solve a system of time dependent differential problems. If a system
     *  of time independent problem is to be solved through a splitting method,
     *  use the class \c dcp::EquationSystem (you need to write the ufl files 
     *  representing the splitting anyway).
     *  It is an abstract class, it only provides the basic interface to 
     *  all splitting methods.
     */ 
    class AbstractSplittingMethod
    {
        // ---------------------------------------------------------------------------------------------//

        public:
            /************************* CONSTRUCTORS ********************/
            //!  Constructor 
            /*!
             *  \param functionSpaces the function spaces of the various differential problems that will be stored in 
             *  the protected member \c differentialSystem_
             */
            AbstractSplittingMethod (const std::vector<std::shared_ptr <dolfin::FunctionSpace>> functionSpaces);

            /************************* DESTRUCTOR ********************/
            //! Destructor
            /*! Default destructor, since members of the class are trivially destructible. 
             *  It is declared as pure virtual to make the class abstract. Indeed, one cannot use this class because
             *  the protected member \c differentialSystem_ is not set (or better, it is set to \c nullptr and thus
             *  using the class methods would result in a segmentation fault). On the other hand, all methods can
             *  be defined here and inherited in derived classes. That's why we defined a pure virtual destructor.
             *  Note that it still needs to be defined, as it will be called anyway when the base class is destroyed.
             */
            virtual ~AbstractSplittingMethod () = 0;


            /********************** GETTERS ***********************/
            //! Get the system stored in the protected member \c differentialSystem_ [1], const version
            /*! 
             *  \return a pointer to the system
             */
            virtual const dcp::TimeDependentEquationSystem& system () const;

            //! Get the system stored in the protected member \c differentialSystem_ [2], non const version
            /*! 
             *  \return a pointer to the system
             */
            virtual dcp::TimeDependentEquationSystem& system ();

            //! Access problem in \c differentialSystem_ with given name [1] (read only)
            /*!
             *  \param name name of the problem to be accessed. 
             *  
             *  \return a reference to the problem
             */
            virtual const dcp::TimeDependentProblem& problem (const std::string& name) const;
            
            //! Access problem in \c differentialSystem_ with given name [2] (read and write)
            /*!
             *  \param name name of the problem to be accessed. 
             *  
             *  \return a reference to the problem
             */
            virtual dcp::TimeDependentProblem& problem (const std::string& name);
            
            
            /********************** SETTERS ***********************/
            //! Set initial solution of the problem with the given name [1].
            //! Wrapper for the method \c setInitialSolution() in \c dcp::TimeDependentProblem
            /*!
             *  \param problemName the name of the problem whose initial solution should be set
             *  \param initialSolution the initial solution itself
             *  \param stepNumber the number of steps to go back to set the initial solution. This is only useful 
             *  for multisteps method, where more than one previous solution is needed to compute the solution at the
             *  current time step. This number, if greater than one, will be concatenated to \c "previous_solution_name"
             *  to form the coefficient name to be set. See the documentation for the constructor of
             *  \c dcp::TimeDependentProblem for more details.
             *  The default value is 1, that means that the input function should be used to set the last solution.
             */
            virtual void setInitialSolution (const std::string& problemName, 
                                             const dolfin::Function& initialSolution,
                                             const unsigned int& stepNumber = 1);
            
            //! Set initial solution of the problem with the given name [1].
            //! Wrapper for the method \c setInitialSolution() in \c dcp::TimeDependentProblem
            /*!
             *  \param problemName the name of the problem whose initial solution should be set
             *  \param initialSolution the initial solution itself
             *  \param stepNumber the number of steps to go back to set the initial solution. This is only useful 
             *  for multisteps method, where more than one previous solution is needed to compute the solution at the
             *  current time step. This number, if greater than one, will be concatenated to \c "previous_solution_name"
             *  to form the coefficient name to be set. See the documentation for the constructor of
             *  \c dcp::TimeDependentProblem for more details.
             *  The default value is 1, that means that the input function should be used to set the last solution.
             */
            virtual void setInitialSolution (const std::string& problemName, 
                                             const dolfin::Expression& initialSolution,
                                             const unsigned int& stepNumber = 1);
            
            //! Set problem coefficients [1].
            //! Wrapper for the method \c setCoefficient() in \c dcp::TimeDependentProblem
            /*!
             *  Parameters are:
             *  \param problemName the name of the problem in which to set the parameter
             *  \param coefficientType used to disambiguate between different member variables to choose which 
             *  coefficient to set
             *  \param coefficientValue value of the coefficient
             *  \param coefficientName string identifying the coefficient to set
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual void setCoefficient (const std::string& problemName,
                                         const std::string& coefficientType, 
                                         const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                         const std::string& coefficientName);

            //! Set problem coefficients [2].
            //! Wrapper for the method \c setCoefficient() in \c dcp::TimeDependentProblem
            /*!
             *  Parameters are:
             *  \param problemName the name of the problem in which to set the parameter
             *  \param coefficientType used to disambiguate between different member variables to choose which
             *  coefficient to set
             *  \param coefficientValue value of the coefficient
             *  \param coefficientNumber integer identifying the coefficient to set
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual void setCoefficient (const std::string& problemName,
                                         const std::string& coefficientType,
                                         const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                         const std::size_t& coefficientNumber);

            //! Set integration subdomains for the forms.
            //! Wrapper for the method \c setIntegrationSubdomain() in \c dcp::TimeDependentProblem
            /*! 
             *  Input arguments are:
             *  \param problemName the name of the problem in which to set the integration subdomain
             *  \param formType used to disambiguate between different member variables to choose which integration
             *  subdomain to set
             *  \param meshFunction the mesh function used to set the integration subdomains
             *  \param subdomainType the type of the subdomains, chosen among those provided by the enumeration
             *  class \c dcp::SubdomainType
             * 
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual void setIntegrationSubdomain (const std::string& problemName,
                                                   const std::string& formType,
                                                   std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                                   const dcp::SubdomainType& subdomainType);

            //! Add Dirichlet boundary condition to the problem [1].
            //! Wrapper for the method \c addDirichletBC() in \c dcp::TimeDependentProblem
            /*!
             *  \param problemName the name of the problem in which to set the Dirichlet boundary condition
             *  \param condition the boundary condition to enforce
             *  \param boundary the boundary on which to enforce the condition
             *  \param bcName the name identifying the boundary condition. If empty,
             *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool addDirichletBC (const std::string& problemName,
                                         const dolfin::GenericFunction& condition, 
                                         const dolfin::SubDomain& boundary,
                                         std::string bcName = "");

            //! Add Dirichlet boundary condition to the problem [2].
            //! Wrapper for the method \c addDirichletBC() in \c dcp::TimeDependentProblem
            /*!
             *  \param problemName the name of the problem in which to set the Dirichlet boundary condition
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
            virtual bool addDirichletBC (const std::string& problemName,
                                         const dolfin::GenericFunction& condition, 
                                         const dolfin::SubDomain& boundary,
                                         const std::size_t& component,
                                         std::string bcName = "");

            //! Add Dirichlet boundary condition to the problem [3].
            //! Wrapper for the method \c addDirichletBC() in \c dcp::TimeDependentProblem
            /*!
             *  \param problemName the name of the problem in which to set the Dirichlet boundary condition
             *  \param condition the boundary condition to enforce
             *  \param boundary the boundary on which to enforce the condition
             *  \param bcName the name identifying the boundary condition. If empty,
             *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool addDirichletBC (const std::string& problemName,
                                         std::shared_ptr<const dolfin::GenericFunction> condition, 
                                         std::shared_ptr<const dolfin::SubDomain> boundary,
                                         std::string bcName = "");

            //! Add Dirichlet boundary condition to the problem [4].
            //! Wrapper for the method \c addDirichletBC() in \c dcp::TimeDependentProblem
            /*!
             *  \param problemName the name of the problem in which to set the Dirichlet boundary condition
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
            virtual bool addDirichletBC (const std::string& problemName,
                                         std::shared_ptr<const dolfin::GenericFunction> condition, 
                                         std::shared_ptr<const dolfin::SubDomain> boundary,
                                         const std::size_t& component,
                                         std::string bcName = "");

            //! Add Dirichlet boundary condition to the problem [5].
            //! Wrapper for the method \c addDirichletBC() in \c dcp::TimeDependentProblem
            /*!
             *  \param problemName the name of the problem in which to set the Dirichlet boundary condition
             *  \param dirichletCondition a const reference to the dirichlet boundary condition to be added to the problem
             *  \param bcName the name identifying the boundary condition. If empty,
             *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool addDirichletBC (const std::string& problemName,
                                         const dolfin::DirichletBC& dirichletCondition, 
                                         std::string bcName = "");

            //! Add Dirichlet boundary condition to the problem [6].
            //! Wrapper for the method \c addDirichletBC() in \c dcp::TimeDependentProblem
            /*!
             *  \param problemName the name of the problem in which to set the Dirichlet boundary condition
             *  \param dirichletCondition a rvalue reference to the dirichlet boundary condition to be added to the problem
             *  \param bcName the name identifying the boundary condition. If empty, 
             *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool addDirichletBC (const std::string& problemName,
                                         dolfin::DirichletBC&& dirichletCondition, 
                                         std::string bcName = "");

            //! Remove Dirichlet boundary condition with given name.
            //! Wrapper for the method \c removeDirichletBC() in \c dcp::TimeDependentProblem
            /*!
             *  \param problemName the name of the problem from which the boundary condition should be removed
             *  \param bcName name of the boundary condition to be removed.
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool removeDirichletBC (const std::string& problemName, const std::string& bcName);
            
            //! Add time dependend Dirichlet boundary condition to the problem [1].
            //! Wrapper for the method \c addTimeDependentDirichletBC() in \c dcp::TimeDependentProblem
            /*! 
             *  \param problemName the name of the problem in which to set the Dirichlet boundary condition
             *  \param condition the callable object defining the boundary condition to enforce
             *  \param boundary the callable object defining the boundary on which to enforce the condition
             *  \param bcName the name identifying the boundary condition. If empty,
             *  "dirichlet_condition_<timeDependentDirichletBCsCounter_>" will be used as default name
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool addTimeDependentDirichletBC (const std::string& problemName,
                                                      const dcp::TimeDependentExpression& condition, 
                                                      const dcp::Subdomain& boundary,
                                                      std::string bcName = "");

            //! Add time dependend Dirichlet boundary condition to the problem [2].
            //! Wrapper for the method \c addTimeDependentDirichletBC() in \c dcp::TimeDependentProblem
            /*! 
             *  \param problemName the name of the problem in which to set the Dirichlet boundary condition
             *  \param condition the callable object defining the boundary condition to enforce
             *  \param boundary the callable object defining the boundary on which to enforce the condition
             *  \param component the function space component on which the boundary condition should be imposed. 
             *  \param bcName the name identifying the boundary condition. If empty,
             *  "dirichlet_condition_<timeDependentDirichletBCsCounter_>" will be used as default name
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool addTimeDependentDirichletBC (const std::string& problemName,
                                                      const dcp::TimeDependentExpression& condition, 
                                                      const dcp::Subdomain& boundary,
                                                      const std::size_t& component,
                                                      std::string bcName = "");

            //! Remove time dependent Dirichlet boundary condition with given name. .
            //! Wrapper for the method \c removeTimeDependentDirichletBC() in \c dcp::TimeDependentProblem
            /*!
             *  \param problemName the name of the problem from which the boundary condition should be removed
             *  \param bcName name of the boundary condition to be removed.
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool removeTimeDependentDirichletBC (const std::string& problemName, const std::string& bcName);
            
            /********************** METHODS ***********************/

            //! Apply method
            /*!
             *  Applies splitting method to solve differential problem.
             *  
             *  \param type the solve type requested. It will be passed to the \c solve method of the protected member
             *  \c differentialSystem_. See \c dcp::AbstractEquationSystem documentation for more details on this 
             *  parameter.
             */
            virtual void apply (const std::string& type = "default");
            
            /********************** VARIABLES ***********************/
            //! the problem parameters
            dolfin::Parameters parameters;
            
            // ---------------------------------------------------------------------------------------------//

        protected:
            //! The problems' finite element spaces
            /*! 
             *  Stored as a \c std::shared_ptr because it may be common to more than 
             *  one problem
             */
            std::vector<std::shared_ptr<dolfin::FunctionSpace> > functionSpaces_;

            //! The system of time dependent differential problems.
            std::shared_ptr<dcp::TimeDependentEquationSystem> differentialSystem_;
            
            // ---------------------------------------------------------------------------------------------//

        private:

    };

}
#endif
