/* 
 *  Copyright (C) 2014, Mattia Tamellini, mattia.tamelllini@gmail.com
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

#ifndef SRC_DIFFERENTIALPROBLEM_ABSTRACTDIFFERENTIALPROBLEM_HPP_INCLUDE_GUARD
#define SRC_DIFFERENTIALPROBLEM_ABSTRACTDIFFERENTIALPROBLEM_HPP_INCLUDE_GUARD

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/fem/DirichletBC.h>
#include <dolfin/parameter/Parameters.h>
#include <boost/shared_ptr.hpp>
#include <DifferentialProblem/SubdomainType.hpp>
#include <map>
#include <string>


namespace DCP
{
    /*! \class AbstractDifferentialProblem AbstractDifferentialProblem.hpp
     *  \brief Abstract base class for differential problems. 
     *         
     *  This class contains the basic interface for a differential problem to be
     *  solved with FEniCS library. It is an abstract class, it only provides the
     *  basic interface to all differential problems
     */ 
    class AbstractDifferentialProblem
    {
        // ---------------------------------------------------------------------------------------------//

        public:
            /************************* CONSTRUCTORS ********************/
            //! Default constructor
            AbstractDifferentialProblem () = delete;

            //!  Constructor with shared pointers
            /*!
             *  \param mesh the problem mesh as a const \c boost::shared_ptr to \c dolfin::Mesh
             *  \param functionSpace the problem finite element space as a const \c boost::shared_ptr to 
             *  \c dolfin::FunctionSpace
             *  The stored mesh's and function space's ownership will be shared between the object and the input argument.
             *  The bilinear and linear form will be created too, calling the constructor which takes the function space
             *  as input.
             */
            AbstractDifferentialProblem (const boost::shared_ptr<dolfin::Mesh> mesh, 
                                         const boost::shared_ptr<dolfin::FunctionSpace> functionSpace);


            //! Constructor with references
            /*!
             *  \param mesh the problem mesh as a const \c dolfin::Mesh&
             *  \param functionSpace the problem finite element space as a const \c dolfin::FunctionSpace&
             *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
             *  initialized using the \c new operator and mesh's and functionSpace's copy constructor
             *  The bilinear and linear form will be created too, calling the constructor which takes the function space
             *  as input.
             */
            AbstractDifferentialProblem (const dolfin::Mesh& mesh, 
                                         const dolfin::FunctionSpace& functionSpace);

            //! Constructor with rvalue references
            /*!
             *  \param mesh the problem mesh as a \c dolfin::Mesh&&
             *  \param functionSpace the problem finite element space as a \c dolfin::FunctionSpace&&
             *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
             *  initialized using the \c new operator and mesh's and functionSpace's move constructor
             *  The bilinear and linear form will be created too, calling the constructor which takes the function space
             *  as input.
             */
            AbstractDifferentialProblem (dolfin::Mesh&& mesh, 
                                         dolfin::FunctionSpace&& functionSpace);


            /************************* DESTRUCTOR ********************/
            //! Destructor
            /*! Default destructor, since members of the class are trivially 
             * destructible.
             */
            virtual ~AbstractDifferentialProblem () {};


            /********************** GETTERS ***********************/
            //! Get problem's mesh
            /*! 
             *  \return a const reference to the problem's mesh
             */
            const dolfin::Mesh& mesh () const;

            //! Get problem's finite element space
            /*! 
             *  \return a const reference to the problem's function space
             */
            const dolfin::FunctionSpace& functionSpace () const;

            //! Get const reference to the problem's dirichlet boundary condition with given name
            /*! 
             *  \param bcName the name identifying the boundary condition
             *  \return a const reference to the problem's dirichletBC identified by \c bcName
             */
            const dolfin::DirichletBC& dirichletBC (const std::string& bcName) const;

            //! Get const reference to the problem's dirichlet boundary conditions map
            /*! 
             *  \return a const reference to the problem's \c dirichletBC map
             */
            const std::map<std::string, dolfin::DirichletBC>& dirichletBCs () const;

            //! Get const reference to the problem's solution
            /*!
             *  \return a const reference to the problem's solution
             */
            const dolfin::Function& solution () const;  

            /********************** SETTERS ***********************/
            //! Set problem coefficients [1]
            /*!
             *  This method is meant to be overridden in derived classes, so that we can have a uniform interface
             *  to set coefficients in all hierarchy.
             *  Parameters are:
             *  \param coefficientType used to disambiguate between different member variables to choose which coefficient
             *  to set
             *  \param coefficientValue value of the coefficient
             *  \param coefficientName string identifying the coefficient to set
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual void setCoefficient (const std::string& coefficientType, 
                                         const boost::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                         const std::string& coefficientName) = 0;

            //! Set problem coefficients [2]
            /*!
             *  This method is meant to be overridden in derived classes, so that we can have a uniform interface
             *  to set coefficients in all hierarchy.
             *  Parameters are:
             *  \param coefficientType used to disambiguate between different member variables to choose which coefficient 
             *  to set
             *  \param coefficientValue value of the coefficient
             *  \param coefficientNumber integer identifying the coefficient to set
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual void setCoefficient (const std::string& coefficientType,
                                         const boost::shared_ptr<const dolfin::GenericFunction> coefficientValue,
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
             *  class \c DCP::SubdomainType
             * 
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual void setIntegrationSubdomains (const std::string& formType,
                                                   boost::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                                   const DCP::SubdomainType& subdomainType) = 0;

            //! Add Dirichlet boundary condition to the problem [1]
            /*!
             *  \param dirichletCondition a const reference to the dirichlet boundary condition to be added to the problem
             *  \param bcName the name identifying the boundary condition. If empty,
             *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool addDirichletBC (const dolfin::DirichletBC& dirichletCondition, std::string bcName = "");

            //! Add Dirichlet boundary condition to the problem [2]
            /*!
             *  \param dirichletCondition a rvalue reference to the dirichlet boundary condition to be added to the problem
             *  \param bcName the name identifying the boundary condition. If empty, 
             *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool addDirichletBC (dolfin::DirichletBC&& dirichletCondition, std::string bcName = "");

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
             * in any concrete instance of the class
             */
            virtual void solve () = 0;
            
            //! Clone method [1]
            /*!
             *  \return a pointer to a \c DCP::AbstractDifferentialProblem containing a copy of the object on 
             *  which it is called. 
             */
            virtual DCP::AbstractDifferentialProblem* clone () const = 0;


            /********************** VARIABLES ***********************/
            //! the problem parameters
            dolfin::Parameters parameters;
            
            // ---------------------------------------------------------------------------------------------//

        protected:
            //! The problem mesh
            /*! 
             *  Stored as a \c boost::shared_ptr because it may be common to more than 
             *  one problem
             */
            boost::shared_ptr<dolfin::Mesh> mesh_;

            //! The problem finite element space
            /*! 
             *  Stored as a \c boost::shared_ptr because it may be common to more than 
             *  one problem
             */
            boost::shared_ptr<dolfin::FunctionSpace> functionSpace_;

            //! The Dirichlet's boundary conditions vector
            std::map<std::string, dolfin::DirichletBC> dirichletBCs_;

            //! Solution of the differential problem
            dolfin::Function solution_;
            
            //! Counter of dirichletBC inserted in the protected member map. It is used to create a unique
            //! name for insertion of dirichlet bcs if the input argument \c bcName to \c addDirichletBC() is
            //! left empty
            int dirichletBCsCounter_;

            // ---------------------------------------------------------------------------------------------//

        private:

    };

}
#endif
