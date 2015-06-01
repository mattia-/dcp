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

#ifndef SRC_DIFFERENTIAL_PROBLEMS_ABSTRACTPROBLEM_H_INCLUDE_GUARD
#define SRC_DIFFERENTIAL_PROBLEMS_ABSTRACTPROBLEM_H_INCLUDE_GUARD

#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/fem/DirichletBC.h>
#include <dolfin/parameter/Parameters.h>
#include <dolfin/plot/plot.h>
#include <dolfin/plot/VTKPlotter.h>
#include <dcp/differential_problems/SubdomainType.h>
#include <map>
#include <string>
#include <utility>


namespace dcp
{
    /*! \class AbstractProblem AbstractProblem.h
     *  \brief Abstract base class for differential problems. 
     *         
     *  This class contains the basic interface for a differential problem to be
     *  solved with FEniCS library. It is an abstract class, it only provides the
     *  basic interface to all differential problems
     */ 
    class AbstractProblem
    {
        // ---------------------------------------------------------------------------------------------//

        public:
            /************************* CONSTRUCTORS ********************/
            //! Default constructor
            AbstractProblem () = delete;

            //!  Constructor with shared pointer
            /*!
             *  \param functionSpace the problem finite element space as a const \c std::shared_ptr to 
             *  \c dolfin::FunctionSpace
             *  The stored function space's ownership will be shared between the object and the input argument.
             *  The constructors also sets the following parameters:
             *      - \c "plot_component" the component of the solution to be plotted (if the solution is vectorial).
             *        A negative value stands for all the components. Default value: -1
             *      - \c "plot_title" the title of the plot. Default value: "Solution"
             *      - \c "clone_method" the type of clone desired. It can be either \c "shallow_clone" or 
             *        \c "deep_clone". The former stores a pointer to the mesh and function space in the cloned 
             *        object, the latter copies the actual objects. Default value: \c "shallow_clone"
             */
            AbstractProblem (const std::shared_ptr<dolfin::FunctionSpace> functionSpace);


            //! Constructor with reference
            /*!
             *  \param functionSpace the problem finite element space as a const \c dolfin::FunctionSpace&
             *  The stored function space's ownership will be unique to the object, since the pointer is 
             *  initialized using the \c new operator and functionSpace's copy constructor.
             *  The constructors also sets the following parameters:
             *      - \c "plot_component" the component of the solution to be plotted (if the solution is vectorial).
             *        A negative value stands for all the components. Default value: -1
             *      - \c "plot_title" the title of the plot. Default value: "Solution"
             *      - \c "clone_method" the type of clone desired. It can be either \c "shallow_clone" or 
             *        \c "deep_clone". The former stores a pointer to the mesh and function space in the cloned 
             *        object, the latter copies the actual objects. Default value: \c "shallow_clone"
             */
            AbstractProblem (const dolfin::FunctionSpace& functionSpace);

            //! Constructor with rvalue reference
            /*!
             *  \param functionSpace the problem finite element space as a \c dolfin::FunctionSpace&&
             *  The stored function space's ownership will be unique to the object, since the pointers are 
             *  initialized using the \c new operator and mesh's and functionSpace's move constructor
             *  The constructors also sets the following parameters:
             *      - \c "plot_component" the component of the solution to be plotted (if the solution is vectorial).
             *        A negative value stands for all the components. Default value: -1
             *      - \c "plot_title" the title of the plot. Default value: "Solution"
             *      - \c "clone_method" the type of clone desired. It can be either \c "shallow_clone" or 
             *        \c "deep_clone". The former stores a pointer to the mesh and function space in the cloned 
             *        object, the latter copies the actual objects. Default value: \c "shallow_clone"
             */
            AbstractProblem (dolfin::FunctionSpace&& functionSpace);


            /************************* DESTRUCTOR ********************/
            //! Destructor
            /*! Default destructor, since members of the class are trivially 
             * destructible.
             */
            virtual ~AbstractProblem () {};


            /********************** GETTERS ***********************/
            //! Get problem's mesh
            /*! 
             *  \return a const reference to the problem's mesh
             */
            virtual std::shared_ptr<const dolfin::Mesh> mesh () const;

            //! Get problem's finite element space
            /*! 
             *  \return a const reference to the problem's function space
             */
            virtual std::shared_ptr<dolfin::FunctionSpace> functionSpace () const;

            //! Get const reference to the problem's dirichlet boundary condition with given name
            /*! 
             *  \param bcName the name identifying the boundary condition
             *  \return a const reference to the problem's dirichletBC identified by \c bcName
             */
            virtual const dolfin::DirichletBC& dirichletBC (const std::string& bcName) const;

            //! Get const reference to the problem's dirichlet boundary conditions map
            /*! 
             *  \return a const reference to the problem's \c dirichletBC map
             */
            virtual const std::map<std::string, dolfin::DirichletBC>& dirichletBCs () const;

            //! Get const reference to the problem's solution
            /*!
             *  \param solutionType the type of the solution to return. This allows to return the solution stored
             *  in \c solution_ or the stashed solution stored in \c stashedSolution_ according to the input string.
             *  Possible values:
             *  \li \c "default" the real solution, stored in \c solution_
             *  \li \c "stashed" the stashed solution
             *  
             *  Default value: \c "default"
             *  
             *  \return a const reference to the second field of the last element of the protected member \c solution
             */
            virtual const dolfin::Function& solution (const std::string& solutionType = "default") const;  

            /********************** SETTERS ***********************/
            //! Set problem coefficients [1]
            /*!
             *  This method is meant to be overridden in derived classes, so that we can have a uniform interface
             *  to set coefficients in all hierarchy.
             *  Parameters are:
             *  \param coefficientType used to disambiguate between different member variables to choose which 
             *  coefficient to set
             *  \param coefficientValue value of the coefficient
             *  \param coefficientName string identifying the coefficient to set
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual void setCoefficient (const std::string& coefficientType, 
                                         const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                         const std::string& coefficientName) = 0;

            //! Set problem coefficients [2]
            /*!
             *  This method is meant to be overridden in derived classes, so that we can have a uniform interface
             *  to set coefficients in all hierarchy.
             *  Parameters are:
             *  \param coefficientType used to disambiguate between different member variables to choose which
             *  coefficient to set
             *  \param coefficientValue value of the coefficient
             *  \param coefficientNumber integer identifying the coefficient to set
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual void setCoefficient (const std::string& coefficientType,
                                         const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                         const std::size_t& coefficientNumber) = 0;

            //! Set integration subdomains for the forms
            /*! 
             *  This method is meant to be overridden in derived classes, so that we can have a uniform interface
             *  to set integration subdomains in all hierarchy.
             *  Input arguments are:
             *  \param formType used to disambiguate between different member variables to choose which integration
             *  subdomain to set
             *  \param meshFunction the mesh function used to set the integration subdomains
             *  \param subdomainType the type of the subdomains, chosen among those provided by the enumeration
             *  class \c dcp::SubdomainType
             * 
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual void setIntegrationSubdomain (const std::string& formType,
                                                   std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                                   const dcp::SubdomainType& subdomainType) = 0;

            //! Add Dirichlet boundary condition to the problem [1]
            /*!
             *  \param condition the boundary condition to enforce
             *  \param boundary the boundary on which to enforce the condition
             *  \param bcName the name identifying the boundary condition. If empty,
             *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool addDirichletBC (const dolfin::GenericFunction& condition, 
                                         const dolfin::SubDomain& boundary,
                                         std::string bcName = "");

            //! Add Dirichlet boundary condition to the problem [2]
            /*!
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
            virtual bool addDirichletBC (const dolfin::GenericFunction& condition, 
                                         const dolfin::SubDomain& boundary,
                                         const std::size_t& component,
                                         std::string bcName = "");

            //! Add Dirichlet boundary condition to the problem [3]
            /*!
             *  \param condition the boundary condition to enforce
             *  \param boundary the boundary on which to enforce the condition
             *  \param bcName the name identifying the boundary condition. If empty,
             *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool addDirichletBC (std::shared_ptr<const dolfin::GenericFunction> condition, 
                                         std::shared_ptr<const dolfin::SubDomain> boundary,
                                         std::string bcName = "");

            //! Add Dirichlet boundary condition to the problem [4]
            /*!
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
                                         std::string bcName = "");

            //! Add Dirichlet boundary condition to the problem [5]
            /*!
             *  \param dirichletCondition a const reference to the dirichlet boundary condition to be added to the problem
             *  \param bcName the name identifying the boundary condition. If empty,
             *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool addDirichletBC (const dolfin::DirichletBC& dirichletCondition, 
                                         std::string bcName = "");

            //! Add Dirichlet boundary condition to the problem [6]
            /*!
             *  \param dirichletCondition a rvalue reference to the dirichlet boundary condition to be added to the problem
             *  \param bcName the name identifying the boundary condition. If empty, 
             *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool addDirichletBC (dolfin::DirichletBC&& dirichletCondition, 
                                         std::string bcName = "");

            //! Remove Dirichlet boundary condition with given name
            /*!
             *  \param bcName name of the boundary condition to be removed.
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool removeDirichletBC (const std::string& bcName);
            
            //! This method is meant to be overridden only if needed in the derived classes. It checks for possible
            //! private members to update (e.g. the solver, if the solver method string in the parameters
            //! was changed but the solver itself was not) and performs the updating. In this class, it is implemented
            //! as an empty funciton, so there is no need to override it if the derived class has no need for this
            //! method
            virtual void update ();

            /********************** METHODS ***********************/

            //! Solve method
            /*!
             * Solves differential problem storing the solution in the private member \c solution_.
             * It is a pure virtual method that needs to be overridden
             * in any concrete instance of the class.
             * 
             * \param solveType the solve type requested. It may be useful to differentiate among different behaviours in
             * the derived classes
             */
            virtual void solve (const std::string& solveType = "default") = 0;
            
            //! Copy stashed solution to \c solution_, thus making it the actual solution of the problem
            virtual void applyStashedSolution ();
                
            //! Method to plot the solution
            /*!
             *  \param plotType the type of the plot desired. In this case, it is not very useful, since the only 
             *  possible value for \c plotType is \c "all" (which will simply plot the only solution stored in the 
             *  protected member \c solution_), but it is useful to have the possibility to choose among different
             *  behaviour in derived classes.
             */
            virtual void plotSolution (const std::string& plotType = "all");
            
            //! Clone method
            /*!
             *  \return a pointer to a \c dcp::AbstractProblem containing a copy of the object on 
             *  which it is called. 
             */
            virtual dcp::AbstractProblem* clone () const = 0;


            /********************** VARIABLES ***********************/
            //! the problem parameters
            dolfin::Parameters parameters;
            
            // ---------------------------------------------------------------------------------------------//

        protected:
            //! The problem finite element space
            /*! 
             *  Stored as a \c std::shared_ptr because it may be common to more than 
             *  one problem
             */
            std::shared_ptr<dolfin::FunctionSpace> functionSpace_;

            //! The Dirichlet's boundary conditions. The map associates the bc's name to the bc itself
            std::map<std::string, dolfin::DirichletBC> dirichletBCs_;

            //! Solution of the differential problem. 
            /*! 
             *  It is declared as a vector of pairs so that it may contain the solution and the time at which it was 
             *  computed on different timesteps in the derived class \c dcp::TimeDependentProblem. 
             *  For steady problems, it will just have size 1, and the time (i.e. the first field of the pair) will 
             *  have a placeholder value equal to -1
             *  It is returned through the method \c solution() with input argument \c "default"
             */
            std::vector<std::pair <double, dolfin::Function>> solution_;
            
            //! The stashed solution
            /*! 
             *  Used when one wants to compute a solution but not store it as the actual solution, since it may be 
             *  just a tentative one (for example if one is subiterating). 
             *  It is returned through the method \c solution() with input argument \c "stashed"
             */
            dolfin::Function stashedSolution_;

            //! Counter of dirichletBC inserted in the protected member map. 
            /*! 
             *  It is used to create a unique name for insertion of dirichlet bcs if the input argument \c bcName 
             *  to \c addDirichletBC() is left empty
             */
            int dirichletBCsCounter_;
            
            //! The plotter for the solution of the problem
            std::shared_ptr<dolfin::VTKPlotter> solutionPlotter_;

            // ---------------------------------------------------------------------------------------------//

        private:

    };

}
#endif
