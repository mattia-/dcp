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

#ifndef SRC_PROBLEMS_TIMEDEPENDENTPROBLEM_H_INCLUDE_GUARD
#define SRC_PROBLEMS_TIMEDEPENDENTPROBLEM_H_INCLUDE_GUARD

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
#include <tuple>
#include <dcp/problems/GenericProblem.h>
#include <dcp/factories/LinearSolverFactory.h>
#include <dcp/functions/TimeDependentFunction.h>
#include <dcp/problems/SubdomainType.h>
#include <dcp/expressions/TimeDependentExpression.h>
#include <dcp/subdomains/Subdomain.h>
#include <dcp/time/Time_.h>

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
     *  It inherits publicly from \c GenericProblem
     *  and it extends its functionalities to a concrete differential
     *  problem. The problem to be solved on each timestep is stored as a
     *  <tt> shared_ptr <dcp::GenericProblem> </tt>
     */

    class TimeDependentProblem : public dcp::GenericProblem
    {
        // ---------------------------------------------------------------------------------------------//

        public:
            /******************* TYPEDEFS **********************/
            typedef std::pair <std::string, std::string>
                    TimeDependentCoefficientKey;
            typedef std::shared_ptr <dcp::TimeDependentExpression> // not const because I may have to set its time
                    TimeDependentExpressionCoefficientValue;
            typedef std::shared_ptr <const dcp::TimeDependentFunction>
                    TimeDependentFunctionCoefficientValue;
            typedef std::string
                    TimeDependentDirichletBCKey;
            typedef std::tuple <std::shared_ptr <dcp::TimeDependentExpression>, // not const because I may have to set its time
                                std::shared_ptr <const dcp::Subdomain>,
                                int>
                    TimeDependentExpressionDirichletBCValue;
            typedef std::tuple <std::shared_ptr <const dcp::TimeDependentFunction>,
                                std::shared_ptr <const dcp::Subdomain>,
                                int>
                    TimeDependentFunctionDirichletBCValue;
            typedef std::tuple<double, std::size_t, std::vector<dolfin::Function> >
                    TimeDependentProblemState;


            /******************* CONSTRUCTORS *******************/
            //! Default constructor is deleted. The class is not default constructable.
            TimeDependentProblem () = delete;

            //!  Constructor [1]
            /*!
             *  \param timeSteppingProblem the problem to be solved on each time step
             *  \param time a pointer to the time object that holds the current time value. Note that, since the class
             *  will store a new shared pointer created from the input argument, modifications to the time value made
             *  by any object that share this time pointer (be it another \c TimeDependentProblem or a
             *  \c TimeDependentExpression) are seen by all other objects created using the same pointer. Note also that
             *  the time value stored in the object pointed by \c time will be set to \c startTime (see next
             *  input argument)
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
             *  \param nTimeSchemeSteps the number of time steps involved in the time stepping problem solution.
             *  For example, implicit Euler is a one-step scheme, so \c nTimeSchemeSteps should be set to 1. BDF2,
             *  on the other hand, is a two-steps time scheme, so \c nTimeSchemeSteps should be set to 2.
             *  When the solution vector is created in the constructor it has size equal to \c nTimeSchemeSteps, and the
             *  times stored in the pairs (time, solution) for the initial solutions are
             *  <tt>startTime - (nTimeSchemeSteps - 1) * dt, startTime - (nTimeSchemeSteps - 2) * dt, ..., startTime</tt>.
             *  That is, it is assumed that the initial solutions are associated with times lower than or equal to
             *  \c startTime.
             *  The default value is 1.
             *
             *  The constructors also sets the following parameters:
             *      - \c "problem_type" a string describing the problem. Default value: \c "time_dependent"
             *      - \c "dt_name" the name of the variable representing the time step in the ufl file describing the
             *        problem. Default value: "dt"
             *      - \c "previous_solution_name" the name of the variable representing the solution at the previous
             *        time step in the ufl file describing the problem. Note that if the problem is a multi-step problem
             *        (that is, if \c nTimeSchemeSteps is greater than 1) the value of \c "previous_solution_name" will
             *        be used as a basename and it will be concatenated with a number ranging from 2 to
             *        \c nTimeSchemeSteps.
             *        This means that in a one-step method the solution at the previous step in the ufl file must have
             *        the name stored in \c "previous_solution_name", while in a n-step method the previous solutions
             *        used must be called \c <"previous_solution_name">, \c <"previous_solution_name">_2,
             *        \c <"previous_solution_name">_3 and so on.
             *        Default value for \c "previous_solution_name": "u_old"
             *      - \c "write_interval" controls when the computed solution is stored to file. Basically, the solution
             *        will be saved every \c write_interval time steps. A value less than or equal to 0 means that the
             *        solution should never be saved to file.
             *        Default value: 0
             *      - \c "plot_interval" controls when the computed solution is plotted. Basically, the solution will be
             *        plotted every \c plot_interval time steps. A value less than or equal to 0 means that the
             *        solution should never be plotted.
             *        Default value: 0
             *      - \c "purge_interval" controls when the computed solution is purged. Basically, the solution vector
             *        will only ever contain no more than \c purge_interval elements (but in any case no less than the
             *        number of elements needed to compute the solution at the nex time step, aka \c nTimeSchemeSteps_),
             *        unless a value lower than \c nTimeSchemeSteps_ is used. In that case the solutions vector is never
             *        purged. Notice that this means that the front of the \c std::vector containing the solutions will
             *        be erased (and all the following elements in the vector shifted accordingly), so this is a
             *        potentially inefficient and costly operation if the value used for this parameter is high.
             *        Default value: \c nTimeSchemeSteps_
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
             *  Furthermore, the constructor modifies the \c GenericProblem parameter \c plot_title and sets its
             *  default value to the empty string, so that by default the plot title contains only the value of the
             *  current time when the plot method is called.
             */
            TimeDependentProblem (const std::shared_ptr<dcp::GenericProblem> timeSteppingProblem,
                                  const std::shared_ptr<dcp::Time> time,
                                  const double& startTime,
                                  const double& dt,
                                  const double& endTime,
                                  std::initializer_list<std::string> dtCoefficientTypes,
                                  std::initializer_list<std::string> previousSolutionCoefficientTypes,
                                  const unsigned int& nTimeSchemeSteps = 1);



            //!  Constructor [2]
            /*!
             *  In this constructor, no time object is given, so a new one is created (which will then be
             *  unique to this object, at least as long as another object is created using this object's time)
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
             *  \param nTimeSchemeSteps the number of time steps involved in the time stepping problem solution.
             *  For example, implicit Euler is a one-step scheme, so \c nTimeSchemeSteps should be set to 1. BDF2,
             *  on the other *  hand, is a two-steps time scheme, so \c nTimeSchemeSteps should be set to 2.
             *  The default value is 1.
             *
             *  The constructors also sets the following parameters:
             *      - \c "problem_type" a string describing the problem. Default value: \c "time_dependent"
             *      - \c "dt_name" the name of the variable representing the time step in the ufl file describing the
             *        problem. Default value: "dt"
             *      - \c "previous_solution_name" the name of the variable representing the solution at the previous
             *        time step in the ufl file describing the problem. Note that if the problem is a multi-step problem
             *        (that is, if \c nTimeSchemeSteps is greater than 1) the value of \c "previous_solution_name" will
             *        be used as a basename and it will be concatenated with a number ranging from 2 to
             *        \c nTimeSchemeSteps.
             *        This means that in a one-step method the solution at the previous step in the ufl file must have
             *        the name stored in \c "previous_solution_name", while in a n-step method the previous solutions
             *        used must be called \c <"previous_solution_name">, \c <"previous_solution_name">_2,
             *        \c <"previous_solution_name">_3 and so on.
             *        Default value for \c "previous_solution_name": "u_old"
             *      - \c "write_interval" controls when the computed solution is stored to file. Basically, the solution
             *        will be saved every \c write_interval time steps. A value less than or equal to 0 means that the
             *        solution should never be saved to file.
             *        Default value: 0
             *      - \c "plot_interval" controls when the computed solution is plotted. Basically, the solution will be
             *        plotted every \c plot_interval time steps. A value less than or equal to 0 means that the
             *        solution should never be plotted.
             *        Default value: 0
             *      - \c "purge_interval" controls when the computed solution is purged. Basically, the solution vector
             *        will only ever contain no more than \c purge_interval elements (but in any case no less than the
             *        number of elements needed to compute the solution at the nex time step, aka \c nTimeSchemeSteps_),
             *        unless a value lower than \c nTimeSchemeSteps_ is used. In that case the solutions vector is never
             *        purged. Notice that this means that the front of the \c std::vector containing the solutions will
             *        be erased (and all the following elements in the vector shifted accordingly), so this is a
             *        potentially inefficient and costly operation if the value used for this parameter is high.
             *        Default value: \c nTimeSchemeSteps_
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
             *  Furthermore, the constructor modifies the \c GenericProblem parameter \c plot_title and sets its
             *  default value to the empty string, so that by default the plot title contains only the value of the
             *  current time when the plot method is called.
             */
            TimeDependentProblem (const std::shared_ptr<dcp::GenericProblem> timeSteppingProblem,
                                  const double& startTime,
                                  const double& dt,
                                  const double& endTime,
                                  std::initializer_list<std::string> dtCoefficientTypes,
                                  std::initializer_list<std::string> previousSolutionCoefficientTypes,
                                  const unsigned int& nTimeSchemeSteps = 1);



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
            virtual std::shared_ptr<const dolfin::FunctionSpace> functionSpace () const override;

            //! Get const reference to the problem's dirichlet boundary condition with given name.
            //! Note that time dependent Dirichlet BC can be found in this map along with "stationary" BCs.
            //! When a time dependent Dirichlet BC is retrieved through this method, its time value will be the last
            //! set during the solve process.
            /*!
             *  \param bcName the name identifying the boundary condition
             *  \return a const reference to the problem's dirichletBC identified by \c bcName
             */
            virtual const dolfin::DirichletBC& dirichletBC (const std::string& bcName) const override;

            //! Get const reference to the problem's dirichlet boundary conditions map
            //! Note that time dependent Dirichlet BC can be found in this map along with "stationary" BCs.
            //! When a time dependent Dirichlet BC is retrieved through this method, its time value will be the last
            //! set during the solve process.
            /*!
             *  \return a const reference to the problem's \c dirichletBC map
             */
            virtual const std::map<std::string, dolfin::DirichletBC>& dirichletBCs () const override;

            //! Get reference to the problem's solutions vector (const version). Note that this returns also the time
            //! values associated with the solutions
            /*!
             *  \return a const reference to the problem's time-solution pairs vector
             */
            virtual const dcp::TimeDependentFunction& solutionsVector () const;

            //! Get reference to the problem's solutions vector (non const version). Note that this returns also the
            //! time values associated with the solutions
            /*!
             *  \return a const reference to the problem's time-solution pairs vector
             */
            virtual dcp::TimeDependentFunction& solutionsVector ();

            //! Get shared pointer to the current simulation time
            /*!
             *  \return the shared pointer stored in the protected member \c time_. The shared pointer is returned as
             *  non-const since one may want to use this function to set a different value for the time or to
             *  build a new time dependent object sharing the same time object.
             */
            virtual std::shared_ptr<const dcp::Time> time () const;

            //! Get a reference to the simulation start time
            /*!
             *  \return a reference to \c startTime_
             */
            virtual const double& startTime ();

            //! Get a reference to the simulation time step
            /*!
             *  \return a reference to \c dt_
             */
            virtual const double& dt ();

            //! Get a reference to the simulation end time
            /*!
             *  \return a reference to \c endTime_
             */
            virtual const double& endTime ();

            //! Get const reference to the problem's time stepping problem
            /*!
             *  \return a const reference to the problem's time stepping problem
             */
            virtual dcp::GenericProblem& timeSteppingProblem ();


            /******************* SETTERS *******************/
            //! Set initial solution for the time loop using a \c dolfin::Function.
            //! Note that this will not clear the solutions vector. Instead, it will change the value of the last
            //! solutions stored in \c solution_ which will then be used to advance in time.
            /*!
             *  \param initialSolution the function to be used as initial solution
             *  \param stepNumber the number of steps to go back to set the initial solution. This is mostly useful for
             *  multistep method, where more than one initial solution is needed.
             *  The default value is 1, that means that the input function should be used to set the solution at the
             *  last time step before the one that is being currently considered.
             */
            virtual void setInitialSolution (const dolfin::Function& initialSolution,
                                             const unsigned int& stepNumber = 1);

            //! Set initial solution for the time loop using a \c dolfin::Expression.
            //! Note that this will not clear the solutions vector. Instead, it will change the value of the last
            //! solutions stored in \c solution_ which will then be used to advance in time.
            /*!
             *  \param initialSolution the function to be used as initial solution
             *  \param stepNumber the number of steps to go back to set the initial solution. This is mostly useful for
             *  multistep method, where more than one initial solution is needed.
             *  The default value is 1, that means that the input expression should be used to set the solution at the
             *  last time step before the one that is being currently considered.
             */
            virtual void setInitialSolution (const dolfin::Expression& initialSolution,
                                             const unsigned int& stepNumber = 1);

            //! Set coefficient [1]. Override of virtual function in \c GenericProblem.
            //! This function is used to set the coefficients for the protected member \c timeSteppingProblem_.
            //! Note that this function only wraps a call to the \c timeSteppingProblem_ \c setCoefficient function.
            /*!
             *  Possible values for \c coefficientType are:
             *  \li \c bilinear_form to set the coefficient in the bilinear form
             *  \li \c linear_form to set the coefficient in the linear form
             *
             *  See \c GenericProblem documentation for more details on the function.
             */
            virtual void setCoefficient (const std::string& coefficientType,
                                         const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                         const std::string& coefficientName) override;

            //! Set coefficient [2]. Override of virtual function in \c GenericProblem.
            //! This function is used to set the coefficients for the protected member \c timeSteppingProblem_.
            //! Note that this function only wraps a call to the \c timeSteppingProblem_ \c setCoefficient function.
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
            //! This function sets the integration subdomains for the protected member \c timeSteppingProblem_.
            //! Note that this function only wraps a call to the \c timeSteppingProblem_ \c setCoefficient function.
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

            //! Add Dirichlet boundary condition to the problem [1]
            //! Overrides method in \c GenericProblem.
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
            //! Overrides method in \c GenericProblem.
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
            //! Overrides method in \c GenericProblem.
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
            //! Overrides method in \c GenericProblem.
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
            //! Overrides method in \c GenericProblem.
            //! This function sets Dirichlet boundary conditions for the protected member \c timeSteppingProblem_.
            //! Note that this function only wraps a call to the \c timeSteppingProblem_ \c addDirichletBC function.
            /*!
             *  \param dirichletCondition a const reference to the dirichlet boundary condition to be added to the
             *  problem
             *  \param bcName the name identifying the boundary condition. If empty,
             *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
             *
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool addDirichletBC (const dolfin::DirichletBC& dirichletCondition,
                                         std::string bcName = "") override;

            //! Add Dirichlet boundary condition to the problem [6].
            //! Overrides method in \c GenericProblem.
            //! This function sets Dirichlet boundary conditions for the protected member \c timeSteppingProblem_.
            //! Note that this function only wraps a call to the \c timeSteppingProblem_ \c addDirichletBC function.
            /*!
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

            //! Add time dependend Dirichlet boundary condition as expression to the problem [1]
            /*!
             *  This function adds the time dependent expression passed as input arguments to
             *  \c timeDependentExpressionDirichletBCs_ and sets the time dependent Dirichlet boundary condition for the
             *  protected member \c timeSteppingProblem_.
             *  Note that both \c condition and \c boundary will be cloned, so be sure the object passed are indeed
             *  clonable. Expressions and subdomains defined in dcp provide a clone method, but if you decide to derive
             *  those classes to add any functionality be sure to override the clone method.
             *
             *  \param condition the boundary condition to enforce
             *  \param boundary the boundary on which to enforce the condition
             *  \param bcName the name identifying the boundary condition. If empty,
             *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
             *
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool addTimeDependentDirichletBC (const dcp::TimeDependentExpression& condition,
                                                      const dcp::Subdomain& boundary,
                                                      std::string bcName = "");

            //! Add time dependend Dirichlet boundary condition as expression to the problem [2]
            /*!
             *  This function adds the time dependent expression passed as input arguments to
             *  \c timeDependentExpressionDirichletBCs_ and sets the time dependent Dirichlet boundary condition for the
             *  protected member \c timeSteppingProblem_.
             *  Note that both \c condition and \c boundary will be cloned, so be sure the object passed are indeed
             *  clonable. Expressions and subdomains defined in dcp provide a clone method, but if you decide to derive
             *  those classes to add any functionality be sure to override the clone method.
             *
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
            virtual bool addTimeDependentDirichletBC (const dcp::TimeDependentExpression& condition,
                                                      const dcp::Subdomain& boundary,
                                                      const std::size_t& component,
                                                      std::string bcName = "");

            //! Add time dependend Dirichlet boundary condition as expression to the problem [3]
            /*!
             *  This function adds the time dependent expression passed as input arguments to
             *  \c timeDependentExpressionDirichletBCs_ and sets the time dependent Dirichlet boundary condition for the
             *  protected member \c timeSteppingProblem_.
             *  Note that both \c condition and \c boundary will be cloned, so be sure the object passed are indeed
             *  clonable. Expressions and subdomains defined in dcp provide a clone method, but if you decide to derive
             *  those classes to add any functionality be sure to override the clone method.
             *
             *  \param condition the boundary condition to enforce
             *  \param boundary the boundary on which to enforce the condition
             *  \param bcName the name identifying the boundary condition. If empty,
             *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
             *
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool addTimeDependentDirichletBC (const std::shared_ptr<const dcp::TimeDependentExpression> condition,
                                                      const std::shared_ptr<const dcp::Subdomain> boundary,
                                                      std::string bcName = "");

            //! Add time dependend Dirichlet boundary condition as expression to the problem [4]
            /*!
             *  This function adds the time dependent expression passed as input arguments to
             *  \c timeDependentExpressionDirichletBCs_ and sets the time dependent Dirichlet boundary condition for the
             *  protected member \c timeSteppingProblem_.
             *  Note that both \c condition and \c boundary will be cloned, so be sure the object passed are indeed
             *  clonable. Expressions and subdomains defined in dcp provide a clone method, but if you decide to derive
             *  those classes to add any functionality be sure to override the clone method.
             *
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
            virtual bool addTimeDependentDirichletBC (const std::shared_ptr<const dcp::TimeDependentExpression> condition,
                                                      const std::shared_ptr<const dcp::Subdomain> boundary,
                                                      const std::size_t& component,
                                                      std::string bcName = "");


            //! Add time dependend Dirichlet boundary condition as function to the problem [1]
            /*!
             *  This function adds the time dependent function passed as input arguments to
             *  \c timeDependentFunctionDirichletBCs_ and sets the time dependent Dirichlet boundary condition for the
             *  protected member \c timeSteppingProblem_.
             *  Note that \c boundary will be cloned, so be sure the object passed is indeed clonable. Subdomains
             *  defined in dcp provide a clone method, but if you decide to derive that class to add any functionality
             *  be sure to override the clone method. The time dependent function \c condition , on the other hand, will
             *  NOT be cloned; a pointer to the object will be used instead, for efficiency reasons. So make sure the
             *  input object's lifespan is long enough. Note that any change to the condition used will be reflected
             *  on the Dirichlet BC enforced as well.
             *  Also note that the condition is supposed to have at least as many elements as the remaining number of
             *  timestep in the problem when the condition is added, and the time values stored in such elements must
             *  match the times which the problem goes through in the solutions loop
             *
             *
             *  \param condition the boundary condition to enforce
             *  \param boundary the boundary on which to enforce the condition
             *  \param bcName the name identifying the boundary condition. If empty,
             *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
             *
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool addTimeDependentDirichletBC (const dcp::TimeDependentFunction& condition,
                                                      const dcp::Subdomain& boundary,
                                                      std::string bcName = "");

            //! Add time dependend Dirichlet boundary condition as function to the problem [2]
            /*!
             *  This function adds the time dependent function passed as input arguments to
             *  \c timeDependentFunctionDirichletBCs_ and sets the time dependent Dirichlet boundary condition for the
             *  protected member \c timeSteppingProblem_.
             *  Note that \c boundary will be cloned, so be sure the object passed is indeed clonable. Subdomains
             *  defined in dcp provide a clone method, but if you decide to derive that class to add any functionality
             *  be sure to override the clone method. The time dependent function \c condition , on the other hand, will
             *  NOT be cloned; a pointer to the object will be used instead, for efficiency reasons. So make sure the
             *  input object's lifespan is long enough. Note that any change to the condition used will be reflected
             *  on the Dirichlet BC enforced as well.
             *  Also note that the condition is supposed to have at least as many elements as the remaining number of
             *  timestep in the problem when the condition is added, and the time values stored in such elements must
             *  match the times which the problem goes through in the solutions loop
             *
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
            virtual bool addTimeDependentDirichletBC (const dcp::TimeDependentFunction& condition,
                                                      const dcp::Subdomain& boundary,
                                                      const std::size_t& component,
                                                      std::string bcName = "");

            //! Add time dependend Dirichlet boundary condition as function to the problem [3]
            /*!
             *  This function adds the time dependent function passed as input arguments to
             *  \c timeDependentFunctionDirichletBCs_ and sets the time dependent Dirichlet boundary condition for the
             *  protected member \c timeSteppingProblem_.
             *  Note that \c boundary will be cloned, so be sure the object passed is indeed clonable. Subdomains
             *  defined in dcp provide a clone method, but if you decide to derive that class to add any functionality
             *  be sure to override the clone method. The time dependent function \c condition , on the other hand, will
             *  NOT be cloned; a pointer to the object will be used instead, for efficiency reasons. So make sure the
             *  input object's lifespan is long enough. Note that any change to the condition used will be reflected
             *  on the Dirichlet BC enforced as well.
             *  Also note that the condition is supposed to have at least as many elements as the remaining number of
             *  timestep in the problem when the condition is added, and the time values stored in such elements must
             *  match the times which the problem goes through in the solutions loop
             *
             *
             *  \param condition the boundary condition to enforce
             *  \param boundary the boundary on which to enforce the condition
             *  \param bcName the name identifying the boundary condition. If empty,
             *  "dirichlet_condition_<dirichletBCsCounter>" will be used as default name
             *
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool addTimeDependentDirichletBC (const std::shared_ptr<const dcp::TimeDependentFunction> condition,
                                                      const std::shared_ptr<const dcp::Subdomain> boundary,
                                                      std::string bcName = "");

            //! Add time dependend Dirichlet boundary condition as function to the problem [4]
            /*!
             *
             *  This function adds the time dependent function passed as input arguments to
             *  \c timeDependentFunctionDirichletBCs_ and sets the time dependent Dirichlet boundary condition for the
             *  protected member \c timeSteppingProblem_.
             *  Note that \c boundary will be cloned, so be sure the object passed is indeed clonable. Subdomains
             *  defined in dcp provide a clone method, but if you decide to derive that class to add any functionality
             *  be sure to override the clone method. The time dependent function \c condition , on the other hand, will
             *  NOT be cloned; a pointer to the object will be used instead, for efficiency reasons. So make sure the
             *  input object's lifespan is long enough. Note that any change to the condition used will be reflected
             *  on the Dirichlet BC enforced as well.
             *  Also note that the condition is supposed to have at least as many elements as the remaining number of
             *  timestep in the problem when the condition is added, and the time values stored in such elements must
             *  match the times which the problem goes through in the solutions loop
             *
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
            virtual bool addTimeDependentDirichletBC (const std::shared_ptr<const dcp::TimeDependentFunction> condition,
                                                      const std::shared_ptr<const dcp::Subdomain> boundary,
                                                      const std::size_t& component,
                                                      std::string bcName = "");

            //! Remove time dependent Dirichlet boundary condition with given name.
            /*!
             *  This function removes the given Dirichlet boundary conditions from the protected members
             *  \c timeDependentExpressionDirichletBCs_ and \c timeDependentFunctionDirichletBCs_ . It also removes it
             *  from \c timeSteppingProblem_.
             *
             *  \param bcName name of the boundary condition to be removed.
             *
             *  \return boolean flag, with \c true representing success and \c false representing failure
             */
            virtual bool removeTimeDependentDirichletBC (const std::string& bcName);

            //! Add an entry to the protected member map \c timeDependentExpressionCoefficients_ .
            /*
             *  \param coefficientType the type of the coefficient, in a form that will make sense once passed to the
             *  method \c setCoefficient() (that is for example \c "linear_form", \c "bilinear_form" and so on)
             *  \param expression the time dependent expression, whose \c eval() method will be used when setting the
             *  coefficient in \c setTimeDependentCoefficients_().
             *  \param coefficientName the name to identify the coefficient. It will be used in the protected method
             *  \c setTimeDependentCoefficients_(), so it must be the same as the name used in the ufl file
             *
             *  \return \c true if the coefficient was inserted in the map, \c false otherwise
             */
            virtual bool addTimeDependentCoefficient
                (const std::string& coefficientType,
                 const std::shared_ptr<dcp::TimeDependentExpression> expression,
                 const std::string& coefficientName);


            //! Add an entry to the protected member map \c timeDependentFunctionCoefficients_ .
            /*
             *  Note that the coefficient is supposed to have at least as many elements as the remaining number of
             *  timestep in the problem when the condition is added, and the time values stored in such elements must
             *  match the times which the problem goes through in the solutions loop
             *
             *  \param coefficientType the type of the coefficient, in a form that will make sense once passed to the
             *  method \c setCoefficient() (that is for example \c "linear_form", \c "bilinear_form" and so on)
             *  \param function the time dependent function which will be used when setting the coefficient in
             *  \c setTimeDependentCoefficients_().
             *  \param coefficientName the name to identify the coefficient. It will be used in the protected method
             *  \c setTimeDependentCoefficients_(), so it must be the same as the name used in the ufl file
             *
             *  \return \c true if the coefficient was inserted in the map, \c false otherwise
             */
            virtual bool addTimeDependentCoefficient (const std::string& coefficientType,
                                                      const std::shared_ptr<const dcp::TimeDependentFunction> function,
                                                      const std::string& coefficientName);

            //! Remove the selected element from the protected member map \c timeDependentExpressionCoefficients_
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
             *  \return \c true if the time loop has ended, that is if current time is greater than or equal to
             *  \c parameters["end_time"], \c false otherwise
             */
            virtual bool isFinished ();

            //! Clear solutions vector
            /*!
             *  This method clears the protected member \c solution_ and creates zero functions as its first elements,
             *  as if the class had just been created. Protected member \c time_ will be reset to \c startTime_.
             *  Finally, the time value in time dependent Dirichlet BCs will be reset to \c startTime_ as well.
             */
            virtual void clear ();

            //! Save problem's state
            /*!
             *  Save current state of the problem (aka the last \c nTimeSchemeSteps_ solutions, the current time and the
             *  current timestep) in the states map (aka \c states_ ) with the given name.
             *  Solutions are saved in the state's vector in the same order, i.e. the last solution computed (which is
             *  the solution computed at the last timestep) is the last element of the state's vector, the last but one
             *  solution is the last but one element and so on.
             *
             *  \param stateName the name to be used to store the state
             *
             *  \return \c true if the state was saved, \c false otherwise
             */
            virtual bool saveState (const std::string& stateName);

            //! Restore problem's state
            /*!
             *  Restore the problem's state to the one in \c state_ with the given name. This means that the values of
             *  \c time and \c timestep will be reset to those of the state and that the state's vector will be added
             *  at the end of <tt>solution_</tt>, so that the solve process can be restarted.
             *  Note that the elements alreay stored in \c solution_ will not be erased, unless the input argument
             *  \c clear is set to \c true. The state will not be erased from \c states_ .
             *
             *  \param stateName the name of the state to be restored
             *  \param clear boolean switch; if set to \c true , then the entire solutions vector will be cleared before
             *  restoring the state, so that the restored state is the only thing left in \c solutions_ after the
             *  function returns. Default value: false
             *
             *  \return \c true if the state was restored, \c false otherwise
             */
            virtual bool restoreState (const std::string& stateName, const bool& clear = false);

            //! Remove state with given name from saved states
            /*!
             *  \param stateName name of the state to be deleted
             *
             *  \return \c true if state was deleted, \c false otherwise
             */
            virtual bool deleteState (const std::string& stateName);

            //! Reserve memory for solutions vector [1]
            /*!
             *  This method calls \c reserve() on \c solution_ so that enough space is allocated to contain a given
             *  amount of elements without reallocation.
             *  \param nElements number of elements to be passed to \c reserve()
             */
            virtual void reserve (std::size_t nElements);

            //! Reserve memory for solutions vector [2]
            /*!
             *  This method calls \c reserve() on \c solution_ so that enough space is allocated to contain a given
             *  amount of elements without reallocation. In this second version of the method, no input argument is
             *  given, so it automatically selects the best capacity for the solutions vector based on initial and
             *  final time, time step and purge interval
             */
            virtual void reserve ();

            //! Advance time value \c time_.
            /*!
             *  It calls the setter function <tt>time_ -> setT0 ()</tt> using <tt>startTime_ + timestep_ * dt_</tt>
             *  and it increments the protected member \c timestep_
             *
             */
            virtual void advanceTime ();

            //! Solve the problem
            /*!
             *  This method solves the problem defined. It uses the protected members' value to set the problem and then
             *  stores the solution in the private member \c solution_ (which is inherited from the base class),
             *  adding to its existing elements (see method \c clear to delete all elements stored in \c solution_).
             *  The initial value for the time loop will be whatever function is stored as the last element of
             *  \c solution_ when the method is called. To change the initial value, use the method \c setInitialSolution.
             *  Note that the value of the solution on the intermediate time steps will be saved int the protected
             *  member \c solution_.  The solution will be plotted every \c plotInterval time steps.
             *  This value is read from  <tt>parameters ["plot_interval"]</tt>. To change the component of the solution
             *  that will be plotted  set <tt>parameters ["plot_components"]</tt> accordingly. By default, this parameter
             *  is set to "-1" (a placeholder that stands for "all the components").
             *
             *  \param solveType the solution type requested. Possible values are:
             *  \li \c "default" the entire time loop is performed
             *  \li \c "step" only one step of the time loop is performed
             *  \li \c "steady" the problem is solved (once), but time is not incremented
             *  \li \c "stash" the problem is solved (once) and time is not incremented, but solution is stored in
             *  \c stashedSolution_
             *  \li <tt> "clear+default" </tt> \c clear method is called before performing the entire time loop
             *  \li <tt> "clear+step" </tt> \c clear method is called before performing one time step of the time loop
             */
            virtual void solve (const std::string& solveType = "default") override;

            //! Copy stashed solution to \c solution_, thus making it the actual solution of the problem, using the
            //! current time value.
            /*!
             *  Note that if time value is the same as the last stored solution, that solution is erased before
             *  copying the stased solution to \c solution_. Otherwise, it is simply added at the end of the vector
             */
            virtual void applyStashedSolution () override;

            //! Plot the solution.
            /*! Overrides method in \c dcp::GenericProblem to take into account the fact that
             *  \c solution_ is now a vector with size greater than one. It uses the value of the parameter \c pause
             *  to decide whether to stop at each plot or not and the value of the parameter \c plot_title to set
             *  the plot title (\c plot_title will actually be added to the time, which is always plotted in the title).
             *  \param plotType select among different behaviour.
             *  Possible values are:
             *  \li \c "default" : plot all the functions in \c solution_
             *  \li \c "last" : plot only the last solution stored in \c solution_
             *  \li \c "position_<N>" : write to file the <tt>N</tt>-th element stored in \c solution_ . An error is
             *  issued if \c N is greater than \c solution_ 's size
             *  \li \c "stashed" : plot the stashed solution
             */
            virtual void plotSolution (const std::string& plotType = "default") override;

            //! Write the solution to file.
            /*! Overrides method in \c dcp::GenericProblem to take into account the fact that \c solution_ is now a
             *  vector with size greater than one.
             *  \param writeType select among different behaviour.
             *  Possible values are:
             *  \li \c "default" : write to file all the functions in \c solution_
             *  \li \c "last" : write to file only the last solution stored in \c solution_
             *  \li \c "position_<N>" : write to file the <tt>N</tt>-th element stored in \c solution_ . An error is
             *  issued if \c N is greater than \c solution_ 's size
             *  \li \c "stashed" : write to file the stashed solution
             */
            virtual void writeSolutionToFile (const std::string& writeType = "default") override;

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
            virtual dcp::TimeDependentProblem* clone () const override;

            // ---------------------------------------------------------------------------------------------//

        protected:
            /******************** METHODS *********************/
            //! Solve the time stepping problem just once without incrementing \c time_
            /*!
             *  This method solves the time stepping problem one time without incrementing the value stored in \c time_.
             *  It uses the protected members' value to set the problem. The last element of \c solution_ will be used
             *  as the previous time step solution. To change it, use the method \c setInitialSolution.
             */
            virtual void steadySolve_ ();

            //! Solve the time stepping problem just once without incrementing \c time_ and store solution in
            //! \c stashedSolution_
            /*!
             *  This method solves the time stepping problem one time without incrementing the value stored in \c time_.
             *  It uses the protected members' value to set the problem. The last element of \c solution_ will be used
             *  as the previous time step solution. To change it, use the method \c setInitialSolution.
             *  The computed solution will be stored in \c stashedSolution_
             */
            virtual void stashSolve_ ();

            //! Perform one step of the time loop
            /*!
             *  This method performes one step of the time loop. It uses the protected members' value to set the problem
             *  and then stores the solution in the private member \c solution_ (which is inherited from the base
             *  class), adding to its existing elements (see method \c clear to delete any element stored in
             *  \c solution_).  The last element of \c solution_ will be used as the previous time step solution. To
             *  change it, use the method \c setInitialSolution.  The protected member \c time_ will be incremented using
             *  the value stored in \c parameters["dt"]
             */
            virtual void step_ ();

            //! Solve the time problem looping through time from \c startTime_ to \c endTime_
            /*!
             *  This method solves the time problem at every time step between \c startTime_ and \c endTime_ by calling
             *  \c step_() on each time step. It will also save the solution to file and plot it according to the values
             *  stored in \c parameters["write_interval"] and \c parameters["plot_interval"]
             *
             */
            virtual void solveLoop_ ();

            //! Set the time dependent Dirichlet boundary conditions at every step of the solve loop
            /*!
             *  For each \c element in \c timeDependentExpressionDirichletBCs_ , it will set the boundary condition's time
             *  using \c t_ and update the bc stored in timeSteppingProblem_ calling \c removeDirichletBC() and
             *  \c addDirichletBC() on \c timeSteppingProblem_ itself. Note that we have to do this since
             *  \c addDirichletBC() always stores a *COPY* of the \c dolfin::DirichletBC object it takes as input
             *  (and we chose to do so in order to allow the use of temporary objects when adding Dirichlet BCs to
             *  the problem), so we cannot store a reference to the \c timeDependentExpressionDirichletBCs_ object in
             *  \c dirichletBCs_ and just change its time value. Also, this would prevent the ability to perform
             *  other operations in on \c timeSteppingProblem_ when modifying its Dirichlet BCs (and this is bad for
             *  example if \c timeSteppingProblem_ is a \c dcp::LinearProblem, since we need to reassemble the system
             *  in this case - see the documentation for \c dcp::LinearProblem , and in particular the role of
             *  the parameter \c "system_is_assembled" )
             */
            virtual void setTimeDependentDirichletBCs_ ();

            //! Reset the time dependent Dirichlet boundary condition pointed by the given iterator [1]
            /*!
             *  This method removes the Dirichlet boundary condition from \c timeSteppingProblem_
             *  and replaces it with a new one with the same name but value updated to the new value of \c time_
             *  \param bcIterator the iterator pointing to the bc that should be replaced
             */
            virtual void resetTimeDependentDirichletBC_
                (std::map <dcp::TimeDependentProblem::TimeDependentDirichletBCKey,
                           dcp::TimeDependentProblem::TimeDependentExpressionDirichletBCValue>::iterator bcIterator);

            //! Reset the time dependent Dirichlet boundary condition pointed by the given iterator [2]
            /*!
             *  This method removes the Dirichlet boundary condition from \c timeSteppingProblem_
             *  and replaces it with a new one with the same name but value updated to the new value of \c time_
             *  \param bcIterator the iterator pointing to the bc that should be replaced
             */
            virtual void resetTimeDependentDirichletBC_
                (std::map <dcp::TimeDependentProblem::TimeDependentDirichletBCKey,
                           dcp::TimeDependentProblem::TimeDependentFunctionDirichletBCValue>::iterator bcIterator);

            //! Set the time dependent coefficients at every step of the solve loop
            /*!
             *  For each \c element in \c timeDependentExpressionCoefficients_ , it will set the coefficient's time
             *  using \c t_ and call \c setCoefficient() ; then for each \c element in
             *  \c timeDependentFunctionCoefficients_ it will call \c setCoefficient()
             */
            virtual void setTimeDependentCoefficients_ ();

            //! Set the previous solution coefficients at every step of the solve loop
            virtual void setPreviousSolutionsCoefficients_ ();

            //! Delete elements from \c solution_ according to \c parameters["purge_interval"]
            virtual void purgeSolutionsVector_ ();

            //! Write given solution to given file. Internal use only
            /*!
             *  \param writer the file writer to be used
             *  \param function the function to be written to file
             *  \param filename name of the file, in case \c writer is \c nullptr
             *  \param t current time value
             */
            virtual void write_ (std::shared_ptr<dolfin::File>& writer,
                                 const std::shared_ptr<const dolfin::Function> function,
                                 const std::string& filename,
                                 const double& t);

            //! Write given solution to given binary file. Internal use only
            /*!
             *  \param writer the binary writer to be used
             *  \param function the function to be written to file
             *  \param filename name of the file, in case \c writer is \c nullptr
             *  \param t current time value
             */
            virtual void write_ (std::shared_ptr<dolfin::HDF5File>& writer,
                                 const std::shared_ptr<const dolfin::Function> function,
                                 const std::string& filename,
                                 const double& t);


            /******************** VARIABLES *********************/
            //! The problem to be solved on each time step
            std::shared_ptr <dcp::GenericProblem> timeSteppingProblem_;

            //! Pointer to the object keeping time for this problem
            std::shared_ptr <dcp::Time> time_;

            //! The start time of the simulation
            double startTime_;

            //! The time step of the simulation
            double dt_;

            //! The end time of the simulation
            double endTime_;

            //! The current time step along the solution loop
            std::size_t timestep_;

            //! A map containing the time dependent expressions needed for the time dependent problem
            /*!
             *  The coefficient is identified by a \c pair containing its name and its type and is associated in the map
             *  to a pointer to \c dcp::TimeDependentExpression
             */
            std::map <dcp::TimeDependentProblem::TimeDependentCoefficientKey,
                      dcp::TimeDependentProblem::TimeDependentExpressionCoefficientValue>
                timeDependentExpressionCoefficients_;

            //! The time dependent Dirichlet's boundary conditions expressed as time dependent expressions
            /*!
             *  The map associates the bc's name to the tuple <time dependent expression, boundary, solution component>
             *  identifying the condition itself. If the condition should be enforced on all the function space
             *  components, \c -1 is used as a placeholder
             */
            std::map <dcp::TimeDependentProblem::TimeDependentDirichletBCKey,
                      dcp::TimeDependentProblem::TimeDependentExpressionDirichletBCValue>
                timeDependentExpressionDirichletBCs_;


            //! A map containing the time dependent functions needed for the time dependent problem
            /*!
             *  The coefficient is identified by a \c pair containing its name and its type and is associated in the map
             *  to a pointer to a \c dcp::TimeDependentFunction (which contains the coefficient itself).
             *  Unlike \c timeDependentExpressionCoefficients_ , in this case the coefficient (which is of type
             *  \c dcp::TimeDependentFunction ) is supposed to have a  different element for each timestep in the
             *  problem (initial and final time included), and the time values stored in such elements must match the
             *  times which the problem goes through in the solutions loop. The right element is selected by using the
             *  protected member \c timestep_ .
             *
             */
            std::map <dcp::TimeDependentProblem::TimeDependentCoefficientKey,
                      dcp::TimeDependentProblem::TimeDependentFunctionCoefficientValue>
                timeDependentFunctionCoefficients_;

            //! The time dependent Dirichlet's boundary conditions expressed as time dependent functions
            /*!
             *  The map associates the bc's name to the tuple <time dependent function, boundary, solution component>
             *  identifying the condition itself.
             *  If the condition should be enforced on all the function space components, \c -1 is used as a
             *  placeholder.
             *  Unlike \c timeDependentExpressionDirichletBCs_ , in this case the condition (which is of type
             *  \c dcp::TimeDependentFunction ) is supposed to have a different element for each timestep in the
             *  problem (initial and final time included), and the time values stored in such elements must match the
             *  times which the problem goes through in the solutions loop. The right element is selected by using the
             *  protected member \c timestep_ .
             */
            std::map <dcp::TimeDependentProblem::TimeDependentDirichletBCKey,
                      dcp::TimeDependentProblem::TimeDependentFunctionDirichletBCValue>
                timeDependentFunctionDirichletBCs_;

            //! The number of steps in the time scheme. For example, implicit Euler method is a 1-step scheme, while
            //! BDF2 is a 2-steps scheme.
            unsigned int nTimeSchemeSteps_;

            //! Counter of time dependent dirichletBC inserted in the protected member map. It is used to create a
            //! unique name for insertion of dirichlet bcs if the input argument \c bcName to
            //! \c addTimeDependentDirichletBC() is left empty
            int timeDependentDirichletBCsCounter_;

            //! States of the problem
            /*!
             *  A map that contains all the saved states of the problem. Each state is identified by a name and contains
             *  a \c double representing the time, a \c std::size_t containig the timestep and a vector of
             *  <tt>dolfin::Function</tt>s representing the solutions at that time. More than one solution may be needed
             *  in case the time advancement uses a multistep scheme.
             *  This should in theory allow the solve process to be restarted from the saved states.
             *  States are saved by calling \c saveState() and restored by calling \c restoreState().
             */
            std::map<std::string, dcp::TimeDependentProblem::TimeDependentProblemState> states_;

            // ---------------------------------------------------------------------------------------------//

        private:
    };
}

#endif
