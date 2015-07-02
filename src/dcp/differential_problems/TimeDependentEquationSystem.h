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

#ifndef SRC_DIFFERENTIAL_PROBLEMS_TIMEDEPENDENTEQUATIONSYSTEM_H_INCLUDE_GUARD
#define SRC_DIFFERENTIAL_PROBLEMS_TIMEDEPENDENTEQUATIONSYSTEM_H_INCLUDE_GUARD

#include <dcp/differential_problems/GenericEquationSystem.h>
#include <dcp/differential_problems/TimeDependentProblem.h>
#include <dcp/time/Time.h>
#include <map>
#include <tuple>
#include <memory>
#include <iostream>
#include <utility>
#include <vector>
#include <functional>

namespace dcp
{
    /*! \class TimeDependentEquationSystem TimeDependentEquationSystem.h
     *  \brief Class for multi-variable and multi-equation coupled time dependent system
     *  
     *  This class derives from GenericEquationSystem and expands its functionalities by defining the
     *  solve method for time dependent systems.
     */
    class TimeDependentEquationSystem : public dcp::GenericEquationSystem
    {
        // ---------------------------------------------------------------------------------------------//  

        public:
            /******************* TYPEDEFS **********************/
            typedef std::tuple <std::string, int, int> PreviousSolutionLinkValue;
            typedef std::pair <LinkKey, PreviousSolutionLinkValue> PreviousSolutionLink;
            

            /******************* CONSTRUCTORS ******************/
            //! Default constructor 
            /*!
             *  \param time the time object for the simulation. It will be compared to the time object stored in the
             *  time dependent problems that are added to the system through the method \c addProblem()
             *  \param startTime the start time for the simulation. It will be compared to the start time stored in the
             *  time dependent problems that are added to the system through the method \c addProblem()
             *  \param dt the time step for the simulation. It will be compared to the time step stored in the
             *  time dependent problems that are added to the system through the method \c addProblem()
             *  \param endTime the end time for the simulation. It will be compared to the end time stored in the
             *  time dependent problems that are added to the system through the method \c addProblem()
             */
            TimeDependentEquationSystem (const std::shared_ptr<dcp::Time>& time,
                                         const double& startTime,
                                         const double& dt,
                                         const double& endTime);
            
            
            /******************* DESTRUCTOR *******************/
            //! Default destructor
            virtual ~TimeDependentEquationSystem () = default;
            

            /******************** METHODS *********************/
            //! Add problem to the map of problems to be solved [1]
            /*!
             *  The parameters are:
             *  \param problemName the problem name
             *  \param problem a const reference to an \c GenericProblem. 
             *  The class will make a copy of the input problem calling the method \c clone().
             *  The problem's name is inserted at the end of \c solveOrder_
             *  This method overrides that in base class, since we need to change the behaviour when a 
             *  \c dcp::TimeDependentProblem is passed to the function. Indeed, if this is the case, the function will
             *  check if the protected members \c time_, \c startTime_, \c endTime_ and \c dt_ are the same as those
             *  stored in the time dependent problem given as input. This means that all time dependent problems 
             *  that one wants to store in an equation system must share the same \c dcp::Time object and have the same
             *  values for start time, end time and time step.
             *  NB: since this method will invoke \c clone() on problem, clone method MUST be \c shallow_clone 
             *  (or an equivalent method that builds the object sharing the same \c dcp::Time object). In fact, the 
             *  system will increment time JUST ONCE, so it is pivotal that all the problems in the system share the
             *  same time object. No check on the clone method is performed in the function because that would limit
             *  the usage of this function to a single clone type, while the user may want to derive a new class from
             *  \c dcp::TimeDependentProblem with a new clone method
             */
            virtual void addProblem (const std::string& problemName, dcp::GenericProblem& problem);
            
            //! Add problem to the map of problems to be solved [2]
            /*!
             *  The parameters are:
             *  \param problemName the problem name
             *  \param problem a shared pointer to a \c dcp::GenericProblem. 
             *  The problem's name is inserted at the end of \c solveOrder_
             *  This method overrides that in base class, since we need to change the behaviour when a 
             *  \c dcp::TimeDependentProblem is passed to the function. Indeed, if this is the case, the function will
             *  check if the protected members \c time_, \c startTime_, \c endTime_ and \c dt_ are the same as those
             *  stored in the time dependent problem given as input. This means that all time dependent problems 
             *  that one wants to store in an equation system must share the same \c dcp::Time object and have the same
             *  values for start time, end time and time step.
             */
            virtual void addProblem (const std::string& problemName, 
                                     const std::shared_ptr<dcp::GenericProblem> problem);
            
            //! Adds link between problems' coefficient and solution at a previous time step [1]
            /*!
             *  This function will add to the stored map \c linksToPreviousSolutions_ a \c std::pair created on the
             *  input arguments and perform the actual linking calling \c linkProblemToPreviousSolution.
             *  \param linkFrom identifies the problem whose parameter (passed as second argument to the function) 
             *  should be linked with the solution of the problem identified by the fourth parameter (\c linkTo)
             *  \param linkedCoefficientName identifies the coefficient to be linked with said solution
             *  \param linkedCoefficientType identifies the type of the coefficient, and will be passed to the function
             *  \c setCoefficient (see \c dcp::GenericProblem documentation for more details)
             *  \param linkTo identifies the problem whose solution is linked to the parameter in the problem
             *  identified by the second and the first arguments respectively. No check is performed on the
             *  existence of such problem
             *  \param nStepsBack the number of steps back we should go to to find the right solution for the link
             *  \param forceRelinking boolean value (default \c false). If the tuple identifying the coefficient already
             *  appears in the protected member variable \c problemsLinks_, it will be relinked using the \c std::pair
             *  passed as first argument if \c forceRelinking is true, and not relinked if it is false (but issuing a
             *  warning in this case)
             */
            virtual void addLinkToPreviousSolution (const std::string& linkFrom, 
                                                    const std::string& linkedCoefficientName,
                                                    const std::string& linkedCoefficientType, 
                                                    const std::string& linkTo,
                                                    const int& nStepsBack,
                                                    const bool& forceRelinking = false);

            //! Adds link between problems' coefficient and solution at a previous time step [2]
            /*!
             *  This function will add to the stored map \c linksToPreviousSolutions_ a \c std::pair created on the
             *  input arguments and perform the actual linking calling \c linkProblemToPreviousSolution.
             *  \param linkFrom identifies the problem whose parameter (passed as second argument to the function) 
             *  should be linked with the solution of the problem identified by the fourth parameter (\c linkTo)
             *  \param linkedCoefficientName identifies the coefficient to be linked with said solution
             *  \param linkedCoefficientType identifies the type of the coefficient, and will be passed to the function
             *  \c setCoefficient (see \c dcp::GenericProblem documentation for more details)
             *  \param linkTo identifies the problem whose solution is linked to the parameter in the problem
             *  identified by the second and the first arguments respectively. No check is performed on the
             *  existence of such problem
             *  \param nStepsBack the number of steps back we should go to to find the right solution for the link
             *  \param linkToComponent identifies the component of the solution of the problem indentified by \c linkTo
             *  that should be used as coefficient in the problem identified by \c linkFrom
             *  \param forceRelinking boolean value (default \c false). If the tuple identifying the coefficient already
             *  appears in the protected member variable \c problemsLinks_, it will be relinked using the \c std::pair
             *  passed as first argument if \c forceRelinking is true, and not relinked if it is false (but issuing a
             *  warning in this case)
             */
            virtual void addLinkToPreviousSolution (const std::string& linkFrom, 
                                                    const std::string& linkedCoefficientName,
                                                    const std::string& linkedCoefficientType, 
                                                    const std::string& linkTo,
                                                    const int& linkToComponent,
                                                    const int& nStepsBack,
                                                    const bool& forceRelinking = false);
            
            //! Remove link between problems' coefficient and solution at a previous time step
            /*!
             *  Removes the link identified by the input arguments from the protected member 
             *  \c linksToPreviousSolutions_ .
             *  The input arguments will be used to create an object of \c dcp::GenericProblem::LinkKey to use
             *  to erase the corresponding entry from \c linksToPreviousSolutions_ .
             *  
             *  \return \c true if the link was removed, \c false otherwise
             */
            bool removeLinkToPreviousSolution (const std::string& linkFrom, 
                                               const std::string& linkedCoefficientName, 
                                               const std::string& linkedCoefficientType);
            
            
            //! Check if system time loop is finished. It basically calls the function \c isFinished() on every problem
            //! stored in \c storedProblems_ and checks if the number of problems whose time loop has ended is equal 
            //! to the size of \c storedProblems_
            /*!
             *  \return true if all problems' time loops have ended, false otherwise
             */
            virtual bool isFinished ();
            
            //! Access problem with given name [1] (read only). 
            //! Overrides method in base class \c dcp::GenericEquationSystem. In this overridden method we return a
            //! reference to \c dcp::TimeDependentProblem since any problem inserted in a time dependent system will
            //! surely be a time dependent problem
            /*!
             *  \param name name of the problem to be accessed. If the name is not found, the function prints an
             *  error message through the function \c dolfin::dolfin_error()
             *  \return a reference to the problem
             */
            virtual const dcp::TimeDependentProblem& operator[] (const std::string& name) const override;
            
            //! Access problem with given name [2] (read and write)
            /*!
            //! Overrides method in base class \c dcp::GenericEquationSystem. In this overridden method we return a
            //! reference to \c dcp::TimeDependentProblem since any problem inserted in a time dependent system will
            //! surely be a time dependent problem
             *  \param name name of the problem to be accessed. If the name is not found, the function prints an
             *  error message through the function \c dolfin::dolfin_error()
             *  \return a reference to the problem
             */
            virtual dcp::TimeDependentProblem& operator[] (const std::string& name) override;
            
            //! Access problem with given position in vector \c solveOrder_ [1] (read only)
            //! Overrides method in base class \c dcp::GenericEquationSystem. In this overridden method we return a
            //! reference to \c dcp::TimeDependentProblem since any problem inserted in a time dependent system will
            //! surely be a time dependent problem
            /*!
             *  \param position position of the problem to be accessed in the private member vector \c solveOrder_. 
             *  If \c position is greater than vector size, the function prints an error message 
             *  through the function \c dolfin::dolfin_error()
             *  \return a reference to the problem
             */
            virtual const dcp::TimeDependentProblem& operator[] (const std::size_t& position) const override;
            
            //! Access problem with given position in vector \c solveOrder_ [2] (read and write)
            //! Overrides method in base class \c dcp::GenericEquationSystem. In this overridden method we return a
            //! reference to \c dcp::TimeDependentProblem since any problem inserted in a time dependent system will
            //! surely be a time dependent problem
            /*!
             *  \param position position of the problem to be accessed in the private member vector \c solveOrder_. 
             *  If \c position is greater than vector size, the function prints an error message 
             *  through the function \c dolfin::dolfin_error()
             *  \return a reference to the problem
             */
            virtual dcp::TimeDependentProblem& operator[] (const std::size_t& position) override;
            
            //! Return the time object of the simulation
            std::shared_ptr<dcp::Time> time () const;
            
            //! Return the start time of the simulation
            const double& startTime () const;
            
            //! Return the time step of the simulation
            const double& dt () const;
            
            //! Return the end time of the simulation
            const double& endTime () const;
             
            //! Advance time value \c time_. It just calls the increment function <tt>time_ -> add ()</tt> 
            //! with <tt>parameters ["dt"]</tt> as input argument.
            //! This allows us to automatically have a backwards time dependent problem if \c dt_ is negative.
            virtual void advanceTime ();
            
            //! Solve all the problems in the order specified by the private member \c solveOrder_. 
            /*! 
             *  The single problems will be solved calling the \c solve method with \c solveType argument equal 
             *  to \c "step"
             */
            virtual void solve () override;
            
            //! Solve the problem corresponding to the name given (once)
            /*!
             *  The problem will be solved calling the \c solve method with \c type argument equal to \c "steady".
             *  \param problemName a string identifying the problem to be solved. If no problem with that name
             *  is found, a warning is issued
             */
            virtual void solve (const std::string& problemName) override;
            
        // ---------------------------------------------------------------------------------------------//  

        protected:
            //! Performs the actual linking to previous solutions
            /*!
             *  Being a protected member, this method is just called from library functions. 
             *  \param link a \c std::pair of the type contained by the protected member \c linksToPreviousSolutions_
             */
            virtual void linkProblemToPreviousSolution (const PreviousSolutionLink& link);
            
            //! The map of the old solutions needed for the system
            /*!
             *  A system may need to have access to old solutions of the single time dependent problems which it is made 
             *  of.
             *  For example one may need to use not only the last solution computed but also the solution at the previous 
             *  time steps, or the solution from two time steps ago. This is what this map is for.
             *  This map associates a <tt> std::tuple<std::string, std::string, std::string></tt> 
             *  and a <tt>std::pair <std::string, int></tt>, where:
             *  \li the first \c string contains the name of the problem whose coefficient should be linked against 
             *  some other problem's solution
             *  \li the second \c string contains the type of such coefficient, in a form that can be passed to
             *  \c dcp::GenericProblem::setCoefficients
             *  \li the third \c string contains the name of said coefficient in the problem
             *  \li the fourth \c string contains the name of the problem whose solution should be used to set the 
             *  coefficient identified by the first three strings
             *  \li the fifth field, which is an \c int, defines which component of the solution of the problem
             *  identified by the fourth \c string should be used in the link. If the whole solution should be used, 
             *  \c -1 is used as a placeholder
             *  \li finally, the last field, which is an \c int, defines the number of timesteps we should go back to to 
             *  find the function we want to link to (that is, if such \c int is equal to 0 the solution at the current
             *  timestep is used, if -1 the solution at the previous timestep, if -2 the soltion two timesteps ago and
             *  so on) 
             */
            std::map <LinkKey, PreviousSolutionLinkValue> linksToPreviousSolutions_;
            
            //! The time of the system
            std::shared_ptr<dcp::Time> time_;
            
            //! Start time for the simulation
            double startTime_;
            
            //! Time step
            double dt_;
            
            //! End time for the simulation
            double endTime_;
    };
}

#endif
