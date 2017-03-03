/* 
 *  Copyright (C) 2015, Ivan Fumagalli, ivan.fumagalli.if@gmail.com
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
#include <dolfin/plot/VTKPlotter.h>
#include <vector>
#include <string>
#include <memory>
#include <initializer_list>
#include <map>
#include <tuple>
#include <dcp/differential_problems/AbstractProblem.h>
#include <dcp/differential_problems/LinearProblem.h>
#include <dcp/differential_problems/TimeDependentProblem.h>
#include <dcp/factories/LinearSolverFactory.h>
#include <dcp/differential_problems/SubdomainType.h>
#include <dcp/expressions/TimeDependentExpression.h>
#include <dcp/subdomains/Subdomain.h>
#include <math.h>

#include <dcp/differential_problems/MeshManager.h> //"MeshManager.h"
#include <dcp/differential_problems/utilities.h> //"utilities.h"

#include "DefaultPostProcessor.h"

#include "GetPot.h"
extern GetPot inputData;

//TODO: namespace->aegir
namespace Ivan
{

    /*! \class MovingTimeDependentProblem MovingTimeDependentProblem.h
     *  \brief Class for time dependent differential problems with moving domain.
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
     *  and \f$ F \left(v\right) : V \rightarrow \mathds{R} \f$ linear form on the same space.
     *  
     *  It inherits publicly from \c TimeDependentProblem and it extends
     *  its functionalities to a moving-domain differential problem.
     */

    class MovingTimeDependentProblem : public dcp::TimeDependentProblem
    {
        // ---------------------------------------------------------------------------------------------//  

        public:
            /******************* TYPEDEFS **********************/
            typedef std::pair <std::string, std::string>            TimeDependentCoefficientKey;
            typedef std::shared_ptr <dcp::TimeDependentExpression>  TimeDependentCoefficientValue;
            typedef std::string                                     TimeDependentDirichletBCKey;
            typedef std::tuple <std::shared_ptr <dcp::TimeDependentExpression>, std::shared_ptr <dcp::Subdomain>, int> 
                    TimeDependentDirichletBCValue;


            /******************* CONSTRUCTORS *******************/
            //! Default constructor is deleted. The class is not default constructable.
            /*!
             * @see dcp::TimeDependentProblem
             */
            MovingTimeDependentProblem () = delete;

            //!  Constructor
            /*!
             *  \param meshManager the mesh manager, owning and moving the (sole) mesh
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
             *  \param wCoefficientTypes an \c initializer_list of strings containing the types of the forms in which 
             *  the parameter whose name is stored in the member variable \c wName should be set. 
             *      //TODO introduce member wName
             *  These strings will be used to call the function \c setCoefficient
             *  on the member variable \c timeSteppingProblem_, so they need to be suitable for that kind of problem.
             *  See the documentation of the function \c setCoefficient in the class \c LinearProblem a
             *  and \c NonlinearProblem for suitable values for this variable.
             *  The values contained in \c wCoefficientTypes will be saved in the member variable \c parameters.
             *  Note that, during the time stepping process, the \c dolfin::Function stored (maybe temporarily) 
             *  in the protected member \c w_ will be used to set the coefficient in the equation whose name is
             *  stored in the parameter \c wName.  If such \c dolfin::Function has more than one
             *  component (i.e. it is a vector function) all of its components will be used by default. To change this
             *  behaviour, one needs to change the value of the parameter \c time_stepping_solution_component.  This
             *  parameter's default value is -1 (which is a placeholder that stands for "use all the solution's
             *  components", but any negative integer will work), but it can be changed to any non-negative integer to
             *  indicate the specific component the time loop should use.
             *      //TODO adapt actual behaviour to documentation, if needed
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
             *      - \c "mesh_displacement_name" the name of the variable representing the mesh displacement 
             *        time step in the ufl file describing the problem. Default value: "w"
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
             *      - \c "mesh_displacement_coefficient_types" (see input arguments documentation)
             *  Furthermore, the constructor modifies the \c AbstractProblem parameter \c plot_title and sets its 
             *  default value to the empty string, so that by default the plot title contains only the value of the 
             *  current time when the plot method is called.
             */
            MovingTimeDependentProblem (const std::shared_ptr<dcp::MeshManager<> > meshManager,
                                        const std::shared_ptr<dcp::AbstractProblem> timeSteppingProblem,
                                        const double& startTime,
                                        const double& dt,
                                        const double& endTime,
                                        std::initializer_list<std::string> dtCoefficientTypes,
                                        std::initializer_list<std::string> previousSolutionCoefficientTypes,
                                        std::initializer_list<std::string> wCoefficientTypes,
                                        int solCompForALEPb = 0,
                                        const unsigned int& nTimeSchemeSteps = 1);


            /******************* DESTRUCTOR *******************/
            //! Destructor
            /*! 
             *  Default destructor, since members of the class are trivially 
             *  destructible.
             */
            virtual ~MovingTimeDependentProblem ();

            
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

            //! Get problem's mesh manager
            /*!
             *  \return a const reference to the problem's mesh manager
             */
            virtual dcp::MeshManager<> & meshManager () const;


            /******************* SETTERS *******************/
            //! Method to set an initial displacement preserving all the mesh marking performed before
            /*!
             *  \param displacement the displacement to set
             */
            virtual void initializeMesh (dolfin::Expression & displacement);

            //! Method to set the post-processor
            virtual void setPostProcessor (Ivan::DefaultPostProcessor * postProcessor);

            
            /******************* METHODS *******************/

            //! Solve the problem
            /*!
             *  This method solves the problem defined. It uses the protected members' value to set the problem and then
             *  stores the solution in the private member \c solution_ (which is inherited from the base class), 
             *  adding to its existing elements (see method \c clear to delete all elements stored in \c solution_).
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
            
            //! Clone method. Overrides method in \c TimeDependentProblem, but performs exactly the same actions
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
            virtual MovingTimeDependentProblem* clone () const override;

            // ---------------------------------------------------------------------------------------------//

        protected:
            
            //! The mesh manager
            std::shared_ptr<dcp::MeshManager<> > meshManager_;
            // TODO preferirei std::shared_ptr<const dcp::MeshManager<> >, ma poi dovrei mettere const tutti i metodi di dcp::MeshManager e usare mutable...

            std::size_t solCompForALEPb_;

            std::shared_ptr<DefaultPostProcessor> postProcessor_;
            
			      //! Save solution, displacement and other output data, used inside the time loop and thus kept protected.
			      //! TODO : introduce this kind of method in dcp::AbstractProblem ?
			      /*!
			       *	\param tmpSolution the current solution
			       *	\param w the mesh displacement
			       *	\param timStep the current time step
			       */
            virtual void saveDataInFile (const dolfin::Function& tmpSolution,
                                         const int& timeStep) const;
          
            //! DEPRECATED version of #saveDatainFile
            virtual void saveDataInFileOld (const dolfin::Function& tmpSolution,
                                         const dolfin::Function& w,
                                         const int& timeStep) const;
    
            //! Update the problem after mesh changing
            /*! After the mesh has changed, the forms, function spaces and functions defined on it have to be updated.
             *  This update is performed via dolfin::adapt
             */
            virtual void adapt ();

            // ---------------------------------------------------------------------------------------------//

        private:


        friend class DefaultPostProcessor;

    };
}
#endif
