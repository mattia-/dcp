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

#ifndef SRC_DIFFERENTIAL_PROBLEMS_TIMEDEPENDENTPROBLEM_H_INCLUDE_GUARD
#define SRC_DIFFERENTIAL_PROBLEMS_TIMEDEPENDENTPROBLEM_H_INCLUDE_GUARD

#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/la/LinearSolver.h>
#include <dolfin/la/GenericLinearSolver.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/Vector.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/parameter/Parameters.h>
#include <dolfin/plot/VTKPlotter.h>
#include <vector>
#include <string>
#include <memory>
#include <initializer_list>
#include <map>
#include <dcp/differential_problems/AbstractProblem.h>
#include <dcp/factories/LinearSolverFactory.h>
#include <dcp/differential_problems/SubdomainType.h>
#include <dcp/expressions/TimeDependentExpression.h>

namespace dcp
{
    /*! \class TimeDependentProblem TimeDependentProblem.h
     *  \brief Class for time dependent differential problems.
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

    class TimeDependentProblem : public dcp::AbstractProblem
    {
        // ---------------------------------------------------------------------------------------------//  

        public:
            /******************* CONSTRUCTORS *******************/
            //! Default constructor is deleted. The class is not default constructable.
            TimeDependentProblem () = delete;

            //!  Constructor
            /*!
             *  \param timeSteppingProblem the problem to be solved on each time step
             *  \param startTime the initial time for the simulation
             *  \param dt the length of the time step to be used
             *  \param endTime the final time for the simulation
             *  \param dtCoefficientTypes an \c initializer_list of strings containing the types of the forms in which 
             *  the parameter whose name is stored in the member variable \c dtName should be set. 
             *  These strings will be used to call the function \c setCoefficient
             *  on the member variable \c timeSteppingProblem_, so they need to be suitable for that kind of problem.
             *  See the documentation of the function \c setCoefficient in the class \c LinearProblem a
             *  and \c NonlinearProblem for suitable values for this variable.
             *  The values contained in \c dtCoefficientTypes will be saved in the member variable \c parameters.
             *  \param previousSolutionCoefficientTypes an \c initializer_list of strings containing the types of the
             *  forms in which the parameter whose name is stored in the member variable\c previousSolutionName should
             *  be set.  These strings will be used to call the function \c setCoefficient on the member variable \c
             *  timeSteppingProblem_, so they need to be suitable for that kind of problem.  See the documentation of
             *  the function \c setCoefficient in the class \c LinearProblem a and \c NonlinearProblem for suitable
             *  values for this variable.
             *  The values contained in \c previousSolutionCoefficientTypes will be saved in the member variable 
             *  \c parameters.
             *  Note that, during the time stepping process, the \c dolfin::Function stored (maybe temporarily) 
             *  in the protected member \c solution_ will be used to set the coefficient in the equation whose name is
             *  stored in the parameter \c previous_solution_name.  If such \c dolfin::Function has more than one
             *  component (i.e. it is a vector function) all of its components will be used by default. To change this
             *  behaviour, one needs to change the value of the parameter \c time_stepping_solution_component.  This
             *  parameter's default value is -1 (which is a placeholder that stands for "use all the solution's
             *  components", but any negative integer will work), but it can be changed to any non-negative integer to
             *  indicate the specific component the time loop should use.
             *  The constructors also sets the following parameters:
             *      - \c "problem_type" a string describing the problem. Default value: \c "time_dependent"
             *      - \c "dt_name" the name of the variable representing the time step in the ufl file describing the
             *        problem. Default value: "dt"
             *      - \c "previous_solution_name" the name of the variable representing the solution at the previous time
             *        step in the ufl file describing the problem. Default value: "u_old"
             *      - \c "store_interval" the interval of time steps to store the solution. Basically, the solution will
             *        be stored in \c solution_ every \c store_interval time steps. A value less than or equal to 0
             *        means that the solution should never be stored. Default value: 1
             *      - \c "plot_interval" the interval of time steps to plot the solution. Basically, the solution will
             *        be plotted every \c plot_interval time steps. A value less than or equal to 0 means that the 
             *        solution should never be plotted. Default value: 0
             *      - \c "time_stepping_solution_component" the component of the solution to be used when advancing the 
             *        time loop (if the solution is vectorial). The function whose name is stored in the parameter
             *        \c "previous_solution_name" will be set using only the component of \c solution_ set in this 
             *        parameter. A negative value stands for all the components. Default value: -1
             *      - \c "pause" if set to \c true, the time stepping loop will stop at each plot and waits for the user
             *        to close the plot window before proceeding. Default value: \c false
             *      - \c "dt_coefficient_types" (see input arguments documentation)
             *      - \c "previous_solution_coefficient_types" (see input arguments documentation)
             *      - \c "previous_solution_is_set_externally" a boolean flag, if set to \c true the problem will not
             *        set the coefficient whose name is stored in the parameter \c "previous_solution_name" at every
             *        time step but will assume that it has already been set externally, for example with a link in
             *        a \c dcp::TimeDependentEquationSystem. Default value: \c false
             *  Furthermore, the constructor modifies the parameter \c plot_title and sets its default value to the
             *  empty string, so that by default the plot title contains only the value of the current time when the 
             *  plot is called.
             */
            TimeDependentProblem (const std::shared_ptr<dcp::AbstractProblem> timeSteppingProblem,
                                  const double& startTime,
                                  const double& dt,
                                  const double& endTime,
                                  std::initializer_list<std::string> dtCoefficientTypes,
                                  std::initializer_list<std::string> previousSolutionCoefficientTypes);


            /******************* DESTRUCTOR *******************/
            //! Destructor
            /*! 
             *  Default destructor, since members of the class are trivially 
             *  destructible.
             */
            virtual ~TimeDependentProblem () {};

            
            /******************* GETTERS *******************/
            //! Get problem's mesh
            /*! 
             *  \return a const reference to the problem's mesh
             */
            virtual std::shared_ptr<const dolfin::Mesh> mesh () const override;

            //! Get problem's finite element space
            /*! 
             *  \return a const reference to the problem's function space
             */
            virtual std::shared_ptr<dolfin::FunctionSpace> functionSpace () const override;

            //! Get const reference to the problem's dirichlet boundary condition with given name
            /*! 
             *  \param bcName the name identifying the boundary condition
             *  \return a const reference to the problem's dirichletBC identified by \c bcName
             */
            virtual const dolfin::DirichletBC& dirichletBC (const std::string& bcName) const override;

            //! Get const reference to the problem's dirichlet boundary conditions map
            /*! 
             *  \return a const reference to the problem's \c dirichletBC map
             */
            virtual const std::map<std::string, dolfin::DirichletBC>& dirichletBCs () const override;

            //! Get const reference to the problem's solution on the last considered time step
            /*!
             *  \return a const reference to the problem's solution on the last considered time step
             */
            virtual const dolfin::Function& solution () const override;  

            //! Get const reference to the problem's solutions vector
            /*!
             *  \return a const reference to the problem's solutions vector
             */
            virtual const std::vector<dolfin::Function>& solutions () const;  
            
            //! Get a vector containing all the solutions stored during the simulations along with the times at which
            //! those solutions were stored
            virtual const std::vector<std::pair <double, std::shared_ptr<const dolfin::Function>> > 
            solutionsWithTimes () const;  
            
            //! Get const reference to the current simulation time
            /*!
             *  \return a const reference to \c t_
             */
            virtual const double& time () const;
            
            //! Get a reference to the simulation start time
            /*!
             *  \return a reference to \c startTime_
             */
            virtual double& startTime ();
            
            //! Get a reference to the simulation time step
            /*!
             *  \return a reference to \c dt_
             */
            virtual double& dt ();
            
            //! Get a reference to the simulation end time
            /*!
             *  \return a reference to \c endTime_
             */
            virtual double& endTime ();
            
            //! Get const reference to the problem's time stepping problem
            /*! 
             *  \return a const reference to the problem's time stepping problem
             */
            virtual dcp::AbstractProblem& timeSteppingProblem ();


            /******************* SETTERS *******************/
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
                                         const std::string& coefficientName) override;

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
                                         const std::size_t& coefficientNumber) override;

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
                                                   const dcp::SubdomainType& subdomainType) override;

            //! Add Dirichlet boundary condition to the problem [1]
            //! Overrides method in \c AbstractProblem.
            //! This function sets Dirichlet boundary conditions for the protected member \c timeSteppingProblem_.
            //! Note that this function only wraps a call to the \c timeSteppingProblem_ \c addDirichletBC function.
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
                                         std::string bcName = "") override;

            //! Add Dirichlet boundary condition to the problem [2]
            //! Overrides method in \c AbstractProblem.
            //! This function sets Dirichlet boundary conditions for the protected member \c timeSteppingProblem_.
            //! Note that this function only wraps a call to the \c timeSteppingProblem_ \c addDirichletBC function.
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
                                         std::string bcName = "") override;

            //! Add Dirichlet boundary condition to the problem [3]
            //! Overrides method in \c AbstractProblem.
            //! This function sets Dirichlet boundary conditions for the protected member \c timeSteppingProblem_.
            //! Note that this function only wraps a call to the \c timeSteppingProblem_ \c addDirichletBC function.
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
                                         std::string bcName = "") override;

            //! Add Dirichlet boundary condition to the problem [4]
            //! Overrides method in \c AbstractProblem.
            //! This function sets Dirichlet boundary conditions for the protected member \c timeSteppingProblem_.
            //! Note that this function only wraps a call to the \c timeSteppingProblem_ \c addDirichletBC function.
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
                                         std::string bcName = "") override;

            //! Add Dirichlet boundary condition to the problem [5]. 
            //! Overrides method in \c AbstractProblem.
            //! This function sets Dirichlet boundary conditions for the protected member \c timeSteppingProblem_.
            //! Note that this function only wraps a call to the \c timeSteppingProblem_ \c addDirichletBC function.
            /*!
             *  \param dirichletCondition a const reference to the dirichlet boundary condition to be added to the problem
             *  \param bcName the name identifying the boundary condition. If empty, 
             *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool addDirichletBC (const dolfin::DirichletBC& dirichletCondition, 
                                         std::string bcName = "") override;

            //! Add Dirichlet boundary condition to the problem [6]. 
            //! Overrides method in \c AbstractProblem.
            //! This function sets Dirichlet boundary conditions for the protected member \c timeSteppingProblem_.
            //! Note that this function only wraps a call to the \c timeSteppingProblem_ \c addDirichletBC function.
            /*!
             *  \param dirichletCondition a rvalue reference to the dirichlet boundary condition to be added to the problem
             *  \param bcName the name identifying the boundary condition. If empty, 
             *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
             *  
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool addDirichletBC (dolfin::DirichletBC&& dirichletCondition, 
                                         std::string bcName = "") override;

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
            virtual bool removeDirichletBC (const std::string& bcName) override;
            
            //! Add an entry to the protected member map \c timeDependentCoefficients_.
            /*
             *  \param coefficientName the name to identify the coefficient. It will be used in the protected method
             *  \c setTimeDependentCoefficients(), so it must be the same as the name used in the ufl file
             *  \param coefficientType the type of the coefficient, in a form that will make sense once passed to the
             *  method \c setCoefficient() (that is for example \c "linear_form", \c "bilinear_form" and so on)
             *  \param expression the time dependent expression, whose \c eval() method will be used when setting the
             *  coefficient in \c setTimeDependentCoefficients()
             *  
             *  \return \c true if the coefficient was inserted in the map, \c false otherwise
             */
            virtual bool addTimeDependentCoefficient (const std::string& coefficientName, 
                                                      const std::string& coefficientType,
                                                      const dcp::TimeDependentExpression& expression);
            
            //! Remove the selected element from the protected member map \c timeDependentCoefficients_
            /*!
             *  \param coefficientName the name of the coefficient to remove
             *  \param coefficientType the type of the coefficient to remove
             *  
             *  \return \c true if the coefficient was removed from the map, \c false otherwise
             */
            virtual bool removeTimeDependentCoefficient (const std::string& coefficientName,
                                                         const std::string& coefficientType);
            
            //! Method to update class members. 
            //! Note that this function only wraps a call to the \c timeSteppingProblem_ \c update function 
            //! (which may be an empty function)
            virtual void update () override;

            
            /******************* METHODS *******************/
            //! Checks if the time loop has ended
            /*!
             *  \return \c true if the time loop has ended, that is if \c t_ >= \c parameters["end_time"],
             *  \c false otherwise
             */
            virtual bool isFinished ();
            
            //! Clear solutions vector
            /*!
             *  This methos clears the protected member \c solution_ and creates a zero function as its first element,
             *  as if the class had just been created. Protected member \c t_ will also be reset to 
             *  <tt>parameters ["start_time"]</tt>
             */  
            virtual void clear ();
             
            //! Perform one step of the time loop
            /*!
             *  This method performes one step of the time loop. It uses the protected members' value to set the problem
             *  and then stores the solution in the private member \c solution_ (which is inherited from the base
             *  class), adding to its existing elements (see method \c clear to delete any element stored in 
             *  \c solution_).  The last element of \c solution_ will be used as the previous time step solution. To
             *  change it, use the method \c setInitialSolution.  The protected member \c t_ will be incremented by 
             *  \c parameters["dt"]
             */
            virtual void step ();

            //! Solve the problem
            /*!
             *  This method solves the problem defined. It uses the protected members' value to set the problem and then
             *  stores the solution in the private member \c solution_ (which is inherited from the base class), 
             *  adding to its existing elements (see method \c clear to delete any element stored in \c solution_).
             *  The initial value for the time loop will be whatever function is stored as the last element of 
             *  \c solution_ when the method is called. To change the initial value, use the method \c setInitialSolution.
             *  Note that the value of the solution on the intermediate time steps will be saved int the protected
             *  member \c solution_. The storing interval is given by <tt>parameters ["store_interval"]</tt>.
             *  The solution will be plotted every \c plotInterval time steps. This value is read from
             *  <tt>parameters ["plot_interval"]</tt>. To change the component of the solution that will be plotted 
             *  set <tt>parameters ["plot_component"]</tt> accordingly. By default, this parameter is set to -1 
             *  (a placeholder that stands for "all the components").
             *  
             *  \param type the solution type requested. Possible values are:
             *  \li \c "default" the entire time loop is performed
             *  \li \c "step" only one step of the time loop is performed
             *  \li <tt> "clear_default" </tt> \c clear method is called before performing the entire time loop
             *  \li <tt> "clear_step" </tt> \c clear method is called before performing one time step of the time loop
             */
            virtual void solve (const std::string& type = "default") override;

            //! Plot method. Overrides the one in \c dcp::AbstractProblem to take into account the fact that 
            //! \c solution_ is now a vector with size greater than one). It uses the value of the parameter \c pause
            //! to decide whether to stop at each plot or not and the value of the parameter \c plot_title to set
            //! the plot title (\c plot_title will actually be added to the time, which is always plotted in the title)
            void plotSolution () override;
            
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
            virtual dcp::TimeDependentProblem* clone () const override;

            // ---------------------------------------------------------------------------------------------//

        protected:
            //! The problem to be solved on each time step
            std::shared_ptr <dcp::AbstractProblem> timeSteppingProblem_;
            
            //! A map containing the time dependent expressions needed for the time dependent problem
            /*!
             *  If for example a problem had an external force that depends on time, its corresponding coefficient in
             *  the time stepping problem would need to be reset at every step, since the time parameter changes. 
             *  The coefficient is identified by a \c pair containing its name and its type and is associated in the map
             *  to a pointer to \c dcp::TimeDependentExpression (the use of the pointer is necessary to call the 
             *  \c eval() function defined in the user-defined class derived from \c dcp::TimeDependentExpression )
             */
            std::map <std::pair <std::string, std::string>, std::shared_ptr <dcp::TimeDependentExpression> > 
                timeDependentCoefficients_;
            
            //! Vector to save the times on which the solution is stored. This way we can return a vector of pairs
            //! <time, solution> in the protected member \c solutionsVector
            std::vector <double> solutionStoringTimes_;
            
            //! The time during the simulation. It takes into account also previous calls to \c solve since it is
            //! not reset after such function is called. It can be reset to the value of the parameter \c t0 
            //! calling the method \c clear
            double t_;
            
            //! The start time of the simulation
            double startTime_;
            
            //! The time step of the simulation
            double dt_;
            
            //! The end time of the simulation
            double endTime_;
            
            //! Method to advance time value \c t_. It just performs the increment <tt>t += parameters ["dt"]</tt>.
            //! This allows us to automatically have a backwards time dependent problem if \c dt is negative.
            /*!
             *  \param previousSolution the solution on the previous time step, to be used when setting the previous 
             *  solution coefficient in the time stepping problem
             */
            virtual void advanceTime (const dolfin::Function& previousSolution);
            
            //! Method to set the time dependent coefficients at every step of the solve loop
            /*
             *  For each \c element in \c timeDependentCoefficients_ , it will set the coefficient time using \c t_ and
             *  calling \c setCoefficient()
             */
            virtual void setTimeDependentCoefficients ();
            
            //! Method to print a warning if \c isFinished() returns \c true. It is just useful to make \c solve()
            //! method clearer to read
            virtual void printFinishedWarning ();
            
            //! Method to store the solution. 
            /*!
             *  \param solution the solution to be stored
             *  \param timeStep the current time step
             *  \param storeInterval the store interval (see constructor documentation)
            */
            void storeSolution (const dolfin::Function& solution, 
                                const int& timeStep, 
                                const int& storeInterval);
            
            //! Method to store the solution on the last time step. 
            /*!
             *  \param solution the solution to be stored
             *  \param timeStep the current time step
             *  \param storeInterval the store interval (see constructor documentation)
            */
            void storeLastStepSolution (const dolfin::Function& solution, 
                                        const int& timeStep, 
                                        const int& storeInterval);
            
            //! Method to plot the solution, used inside the time loop and thus kept protected. It overloads the 
            //! plot method in \c dcp::AbstractProblem, which is still usable (and actually overridden in this class
            //! to take into account the fact that \c solution_ is now a vector with size greater than one)
            /*!
             *  \param solution the solution to be plotted
             *  \param timeStep the current time step
             *  \param plotInterval the plot interval (see constructor documentation)
             *  \param plotComponent the component of the solution to be plotted (see constructor documentation)
             *  \param pause boolean flag, true if the function should wait after plotting
            */
            void plotSolution (dolfin::Function& solution, 
                               const int& timeStep, 
                               const int& plotInterval, 
                               const int& plotComponent,
                               const bool& pause);

            // ---------------------------------------------------------------------------------------------//

        private:
    };
}
#endif
