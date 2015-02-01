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

#ifndef IVAN_MOVING_DIFFERENTIAL_PROBLEMS_TIMEDEPENDENTPROBLEM_H_INCLUDE_GUARD
#define IVAN_MOVING_DIFFERENTIAL_PROBLEMS_TIMEDEPENDENTPROBLEM_H_INCLUDE_GUARD

#include <dolfin/mesh/Mesh.h>
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
#include <differential_problems/AbstractProblem.h>
#include <factories/LinearSolverFactory.h>
#include <differential_problems/SubdomainType.h>
#include "geometry.h"

namespace Ivan
{
    /*! \class MovingTimeDependentProblem MovingTimeDependentProblem.h
     *  \brief Class for linear differential problems.
     *
     *  This class represents problem of the form
     *  \f[
     *      \mbox{Find } u \in V : 
     *      a \left(u \left( t \right), v\right) 
     *      = 
     *      F \left(v\right) 
     *      \ \forall\,v\,\in\,V \ \forall\,t\,\in\,\left[0, T\right]
     *  \f]
     *  with \f$ a \left(u \left(t\right), v\right) : V \times V \rightarrow \mathds{R}\f$ generic form on \f$V\f$
     *  and \f$ L \left(v\right) : V \rightarrow \mathds{R} \f$ linear form on the same space.
     *  
     *  It inherits publicly from \c AbstractProblem
     *  and it extends its functionalities to a concrete differential
     *  problem. The problem to be solved on each timestep is stored as a 
     *  <tt> shared_ptr <dcp::AbstractProblem> </tt>
     */

    class MovingTimeDependentProblem : public dcp::AbstractProblem
    {
        // ---------------------------------------------------------------------------------------------//  

        public:
            /******************* CONSTRUCTORS *******************/
            //! Default constructor is deleted. The class is not default constructable.
            MovingTimeDependentProblem () = delete;

            //!  Constructor with shared pointers
            /*!
             *  \param mesh the problem mesh as a <tt> const std::shared_ptr </tt> to \c dolfin::Mesh
             *  \param functionSpace the problem finite element space as a <tt> const std::shared_ptr </tt> to 
             *  \c dolfin::FunctionSpace
             *  \param startTime the initial time for the simulation
             *  \param dt the length of the time step to be used
             *  \param endTime the final time for the simulation
             *  \param dtCoefficientTypes a vector of strings containing the types of the forms in which 
             *  the parameter whose name is stored in the member variable \c dtName should be set. 
             *  These strings will be used to call the function \c setCoefficient
             *  on the member variable \c timeSteppingProblem_, so they need to be suitable for that kind of problem.
             *  See the documentation of the function \c setCoefficient in the class \c LinearProblem a
             *  and \c NonlinearProblem for suitable values for this variable.
             *  The values contained in \c dtCoefficientTypes will be saved in the member variable \c parameters.
             *  \param previousSolutionCoefficientTypes a vector of strings containing the types of the forms in which 
             *  the parameter whose name is stored in the member variable\c previousSolutionName should be set. 
             *  These strings will be used to call the function \c setCoefficient
             *  on the member variable \c timeSteppingProblem_, so they need to be suitable for that kind of problem.
             *  See the documentation of the function \c setCoefficient in the class \c LinearProblem a
             *  and \c NonlinearProblem for suitable values for this variable.
             *  The values contained in \c previousSolutionCoefficientTypes will be saved in the member variable 
             *  \c parameters.
             *  \param timeSteppingProblem the problem to be solved on each time step. Default value = nullptr
             *  \param storeInterval interval to be used for solution storing. The solution will be stored in the 
             *  private member \c solutions_ every  \c storeInterval time steps. A value less than or equal to 0
             *  means that no intermediate-step solution will be saved. Default value: 1.
             *  \param plotInterval interval to be used for solution plot. The solution will be plotted every 
             *  \c plotInterval time steps. A value less than or equal to 0 means that no plot will pe performed.
             *  Default value: 1.
             *  \param dtName the name of the coefficient representing the time step in the ufl file describing the
             *  problem. Default value: "dt"
             *  \param previousSolutionName the name of the function representing the solution at the previous time step
             *  in the ufl file describing the problem. Default value: "u_old"
             *  Note that, during the time stepping process, the \c dolfin::Function stored in the protected member 
             *  \c solution_ will be used to set the coefficient in the equation whose name is stored in 
             *  \c previousSolutionName. If such \c dolfin::Function has more than one component (i.e. it is a vector
             *  function) all of its component will be used by default. To change this behaviour, one needs to change
             *  the value of the parameter "time_stepping_solution_component", in the public member \c parameters.
             *  This parameter's default value is -1 (which is a placeholder that stands for "use all the solution's
             *  components", but any negative integer will work), but it can be changed to any non-negative 
             *  integer to indicate the specific component the time loop should use.
             */
            MovingTimeDependentProblem 
                (const std::shared_ptr<dolfin::Mesh> mesh, 
                 const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
                 const double& startTime,
                 const double& dt,
                 const double& endTime,
                 const std::vector<std::string>& dtCoefficientTypes,
                 const std::vector<std::string>& previousSolutionCoefficientTypes,
                 const std::shared_ptr<dcp::AbstractProblem> timeSteppingProblem = nullptr,
                 const int& storeInterval = 1,
                 const int& plotInterval = 1,
                 const std::string& dtName = "dt",
                 const std::string& previousSolutionName = "u_old");

            //! Constructor with references
            /*!
             *  \param mesh the problem mesh as a <tt> const dolfin::Mesh& </tt>
             *  \param functionSpace the problem finite element space as a <tt> const dolfin::FunctionSpace& </tt>
             *  \param startTime the initial time for the simulation
             *  \param dt the length of the time step to be used
             *  \param endTime the final time for the simulation
             *  \param dtCoefficientTypes a vector of strings containing the types of the forms in which 
             *  the parameter whose name is stored in the member variable \c dtName should be set. 
             *  These strings will be used to call the function \c setCoefficient
             *  on the member variable \c timeSteppingProblem_, so they need to be suitable for that kind of problem.
             *  See the documentation of the function \c setCoefficient in the class \c LinearProblem a
             *  and \c NonlinearProblem for suitable values for this variable.
             *  The values contained in \c dtCoefficientTypes will be saved in the member variable \c parameters.
             *  \param previousSolutionCoefficientTypes a vector of strings containing the types of the forms in which 
             *  the parameter whose name is stored in the member variable\c previousSolutionName should be set. 
             *  These strings will be used to call the function \c setCoefficient
             *  on the member variable \c timeSteppingProblem_, so they need to be suitable for that kind of problem.
             *  See the documentation of the function \c setCoefficient in the class \c LinearProblem a
             *  and \c NonlinearProblem for suitable values for this variable.
             *  The values contained in \c previousSolutionCoefficientTypes will be saved in the member variable 
             *  \c parameters.
             *  \param timeSteppingProblem the problem to be solved on each time step. Default value = nullptr
             *  \param storeInterval interval to be used for solution storing. The solution will be stored in the 
             *  private member \c solutions_ every  \c storeInterval time steps. A value less than or equal to 0
             *  means that no intermediate-step solution will be saved. Default value: 1.
             *  \param plotInterval interval to be used for solution plot. The solution will be plotted every 
             *  \c plotInterval time steps. A value less than or equal to 0 means that no plot will pe performed.
             *  Default value: 1.
             *  \param dtName the name of the coefficient representing the time step in the ufl file describing the
             *  problem. Default value: "dt"
             *  \param previousSolutionName the name of the function representing the solution at the previous time step
             *  in the ufl file describing the problem. Default value: "u_old"
             *  Note that, during the time stepping process, the \c dolfin::Function stored in the protected member 
             *  \c solution_ will be used to set the coefficient in the equation whose name is stored in 
             *  \c previousSolutionName. If such \c dolfin::Function has more than one component (i.e. it is a vector
             *  function) all of its component will be used by default. To change this behaviour, one needs to change
             *  the value of the parameter "time_stepping_solution_component", in the public member \c parameters.
             *  This parameter's default value is -1 (which is a placeholder that stands for "use all the solution's
             *  components", but any negative integer will work), but it can be changed to any non-negative 
             *  integer to indicate the specific component the time loop should use.
             */
            MovingTimeDependentProblem 
                (const dolfin::Mesh& mesh, 
                 const dolfin::FunctionSpace& functionSpace,
                 const double& startTime,
                 const double& dt,
                 const double& endTime,
                 const std::vector<std::string>& dtCoefficientTypes,
                 const std::vector<std::string>& previousSolutionCoefficientTypes,
                 const std::shared_ptr<dcp::AbstractProblem> timeSteppingProblem = nullptr,
                 const int& storeInterval = 1,
                 const int& plotInterval = 1,
                 const std::string& dtName = "dt",
                 const std::string& previousSolutionName = "u_old");
            

            //! Constructor with rvalue references
            /*!
             *  \param mesh the problem mesh as a <tt> const dolfin::Mesh&& </tt>
             *  \param functionSpace the problem finite element space as a <tt> const dolfin::FunctionSpace&& </tt>
             *  \param startTime the initial time for the simulation
             *  \param dt the length of the time step to be used
             *  \param endTime the final time for the simulation
             *  \param dtCoefficientTypes a vector of strings containing the types of the forms in which 
             *  the parameter whose name is stored in the member variable \c dtName should be set. 
             *  These strings will be used to call the function \c setCoefficient
             *  on the member variable \c timeSteppingProblem_, so they need to be suitable for that kind of problem.
             *  See the documentation of the function \c setCoefficient in the class \c LinearProblem a
             *  and \c NonlinearProblem for suitable values for this variable.
             *  The values contained in \c dtCoefficientTypes will be saved in the member variable \c parameters.
             *  \param previousSolutionCoefficientTypes a vector of strings containing the types of the forms in which 
             *  the parameter whose name is stored in the member variable\c previousSolutionName should be set. 
             *  These strings will be used to call the function \c setCoefficient
             *  on the member variable \c timeSteppingProblem_, so they need to be suitable for that kind of problem.
             *  See the documentation of the function \c setCoefficient in the class \c LinearProblem a
             *  and \c NonlinearProblem for suitable values for this variable.
             *  The values contained in \c previousSolutionCoefficientTypes will be saved in the member variable 
             *  \c parameters.
             *  \param timeSteppingProblem the problem to be solved on each time step. Default value = nullptr
             *  \param storeInterval interval to be used for solution storing. The solution will be stored in the 
             *  private member \c solutions_ every  \c storeInterval time steps. A value less than or equal to 0
             *  means that no intermediate-step solution will be saved. Default value: 1.
             *  \param plotInterval interval to be used for solution plot. The solution will be plotted every 
             *  \c plotInterval time steps. A value less than or equal to 0 means that no plot will pe performed.
             *  Default value: 1.
             *  \param dtName the name of the coefficient representing the time step in the ufl file describing the
             *  problem. Default value: "dt"
             *  \param previousSolutionName the name of the function representing the solution at the previous time step
             *  in the ufl file describing the problem. Default value: "u_old"
             *  Note that, during the time stepping process, the \c dolfin::Function stored in the protected member 
             *  \c solution_ will be used to set the coefficient in the equation whose name is stored in 
             *  \c previousSolutionName. If such \c dolfin::Function has more than one component (i.e. it is a vector
             *  function) all of its component will be used by default. To change this behaviour, one needs to change
             *  the value of the parameter "time_stepping_solution_component", in the public member \c parameters.
             *  This parameter's default value is -1 (which is a placeholder that stands for "use all the solution's
             *  components", but any negative integer will work), but it can be changed to any non-negative 
             *  integer to indicate the specific component the time loop should use.
             */
            MovingTimeDependentProblem 
                (dolfin::Mesh&& mesh, 
                 dolfin::FunctionSpace&& functionSpace,
                 const double& startTime,
                 const double& dt,
                 const double& endTime,
                 const std::vector<std::string>& dtCoefficientTypes,
                 const std::vector<std::string>& previousSolutionCoefficientTypes,
                 const std::shared_ptr<dcp::AbstractProblem> timeSteppingProblem = nullptr,
                 const int& storeInterval = 1,
                 const int& plotInterval = 1,
                 const std::string& dtName = "dt",
                 const std::string& previousSolutionName = "u_old");


            /******************* DESTRUCTOR *******************/
            //! Destructor
            /*! 
             *  Default destructor, since members of the class are trivially 
             *  destructible.
             */
            virtual ~MovingTimeDependentProblem () {};

            
            /******************* GETTERS *******************/
            //! Get const reference to the problem's time stepping problem
            /*! 
             *  \return a const reference to the problem's time stepping problem
             */
            virtual dcp::AbstractProblem& timeSteppingProblem ();

            //! Get const reference to the problem's solution on the last considered time step
            /*!
             *  \return a const reference to the problem's solution on the last considered time step
             */
            virtual const dolfin::Function& solution () const;  

            //! Get const reference to the problem's solutions vector
            /*!
             *  \return a const reference to the problem's solutions vector
             */
            virtual const std::vector<dolfin::Function>& solutionsVector () const;  


            /******************* SETTERS *******************/
            //! Set time stepping problem
            virtual void setTimeSteppingProblem (const std::shared_ptr<dcp::AbstractProblem> timeSteppingProblem);
            
            //! Set initial solution for the time loop using a \c dolfin::Function
            virtual void setInitialSolution (const dolfin::Function& initialSolution);
            
            //! Set initial solution for the time loop using a \c dolfin::Expression
            virtual void setInitialSolution (const dolfin::Expression& initialSolution);
            
            //! Set coefficient [1]. Override of virtual function in \c AbstractProblem.
            //! This function is used to set the coefficients for the protected member \c timeSteppingProblem_.
            //! Note that this function only wraps a call to the \c timeSteppingProblem_ \c setCoefficient function.
            /*!
             *  Possible values for \c coefficientType are:
             *  \li \c bilinear_form to set the coefficient in the bilinear form
             *  \li \c linear_form to set the coefficient in the linear form
             *  
             *  See \c AbstractProblem documentation for more details on the function.
             */
            virtual void setCoefficient (const std::string& coefficientType, 
                                         const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                         const std::string& coefficientName);

            //! Set coefficient [2]. Override of virtual function in \c AbstractProblem.
            //! This function is used to set the coefficients for the protected member \c timeSteppingProblem_.
            //! Note that this function only wraps a call to the \c timeSteppingProblem_ \c setCoefficient function.
            /*!
             *  Possible values for \c coefficientType are:
             *  \li \c bilinear_form to set the coefficient in the bilinear form
             *  \li \c linear_form to set the coefficient in the linear form
             *  
             *  See \c AbstractProblem documentation for more details on the function
             */
            virtual void setCoefficient (const std::string& coefficientType,
                                         const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                         const std::size_t& coefficientNumber);

            //! Set integration subdomains for the forms. Override of virtual function in \c AbstractProblem
            //! This function sets the integration subdomains for the protected member \c timeSteppingProblem_.
            //! Note that this function only wraps a call to the \c timeSteppingProblem_ \c setCoefficient function.
            /*! 
             *  Possible values for \c formType are:
             *  \li \c bilinear_form to set the integration subdomain in the bilinear form
             *  \li \c linear_form to set the integration subdomain in the linear form
             *  
             *  See \c AbstractProblem documentation for more details on the function
             */
            virtual void setIntegrationSubdomains (const std::string& formType,
                                                   std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                                   const dcp::SubdomainType& subdomainType);

            //! Add Dirichlet boundary condition to the problem [1]. Overrides method in \c AbstractProblem
            //! This function sets Dirichlet boundary conditions for the protected member \c timeSteppingProblem_.
            //! Note that this function only wraps a call to the \c timeSteppingProblem_ \c addDirichletBC function.
            /*!
             *  \param dirichletCondition a const reference to the dirichlet boundary condition to be added to the problem
             *  \param bcName the name identifying the boundary condition. If empty, 
             *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool addDirichletBC (const dolfin::DirichletBC& dirichletCondition, std::string bcName = "");

            //! Add Dirichlet boundary condition to the problem [2]. Overrides method in \c AbstractProblem
            //! This function sets Dirichlet boundary conditions for the protected member \c timeSteppingProblem_.
            //! Note that this function only wraps a call to the \c timeSteppingProblem_ \c addDirichletBC function.
            /*!
             *  \param dirichletCondition a rvalue reference to the dirichlet boundary condition to be added to the problem
             *  \param bcName the name identifying the boundary condition. If empty, 
             *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool addDirichletBC (dolfin::DirichletBC&& dirichletCondition, std::string bcName = "");

            //! Remove Dirichlet boundary condition with given position. Overrides method in \c AbstractProblem
            //! This function removes the given Dirichlet boundary conditions from the protected member 
            //! \c timeSteppingProblem_.
            //! Note that this function only wraps a call to the \c timeSteppingProblem_ \c removeDirichletBC function.
            /*!
             *  This method adds to the base class method the setting of parameter \c system_is_assembled to \c false.
             *  \param bcName name of the boundary condition to be removed.
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool removeDirichletBC (const std::string& bcName);
            
            //! Method to update class members. 
            //! Note that this function only wraps a call to the \c timeSteppingProblem_ \c update function 
            //! (which may be an empty function)
            virtual void update ();

            
            /******************* METHODS *******************/
            //! Solve problem
            /*!
             *  This method solves the problem defined. It uses the protected members' value to set the problem and then
             *  stores the solution in the private member \c solution_ (which is inherited from the base class). 
             *  The initial value for the time loop will be whatever function is stored inside of \c solution_ when
             *  the method is called. To change the initial value, use the method \c setInitialSolution.
             *  Note that the value of the solution on the intermediate time steps will be saved int the protected
             *  member \c solutions_
             */
            virtual void solve ();

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
            virtual MovingTimeDependentProblem* clone () const;

            // ---------------------------------------------------------------------------------------------//

        protected:
            //! The problem to be solved on each time step
            std::shared_ptr <dcp::AbstractProblem> timeSteppingProblem_;
            
            //! A vector containing the solutions on the different timesteps, saved every \c storeInterval 
            //! time steps. Note that the value of the variable \c storeInterval is saved in the 
            //! public member \c parameters
            std::vector <dolfin::Function> solutions_;

            //! The mesh manager
            geometry::MeshManager<dolfin::ALE,dolfin::FunctionSpace> meshManager_;

            //! The displacement function
            // TODO : std::shared_ptr<dolfin::GenericFunction> displacement_;
            std::shared_ptr<geometry::MapTgamma> displacement_;
            std::shared_ptr<dolfin::MeshDisplacement> actualDisplacement_;

            // ---------------------------------------------------------------------------------------------//

        private:
    };
}
#endif
