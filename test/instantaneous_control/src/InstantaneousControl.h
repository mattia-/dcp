/* 
 *  Copyright (C) 2017, Ivan Fumagalli, ivan.fumagalli@polimi.it
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

#ifndef SRC_DIFFERENTIAL_PROBLEMS_INSTANTANEOUSCONTROL_H_INCLUDE_GUARD
#define SRC_DIFFERENTIAL_PROBLEMS_INSTANTANEOUSCONTROL_H_INCLUDE_GUARD

//#include <dcp/differential_problems/AbstractEquationSystem.h>
#include <dcp/differential_problems/EquationSystem.h>
#include <dcp/differential_problems/MeshManager.h>
#include <map>
#include <tuple>
#include <memory>
#include <iostream>
#include <utility>
#include <vector>
#include <functional>

namespace dcp
{

class RescaleThroughTime : public dolfin::Expression
{
  public:

    void eval(dolfin::Array<double> &values, const dolfin::Array<double> &x) const
    {
      values[0] = value();
    }

    void setFactor (double scaleFactor)
    {
      scaleFactor_ = scaleFactor;
    }

    /*void setTime (double t)
    {
      t_ = t;
    }*/

    void setN (double n)
    {
      n_ = n;
    }

    void increaseN ()
    {
      ++n_;
    }

    double value() const
    {
      //return scaleFactor_*exp(n_/100.0);
      //return scaleFactor_*sqrt(sqrt(n_/100.0));
      return scaleFactor_;//*n_/100.0 *n_/100.0 *n_/100.0 *n_/100.0;
      //return scaleFactor_*t_;//*t_;
    }

  private:

    double scaleFactor_, n_;
    //double scaleFactor_, t_;
};

    /*! \class InstantaneousControl InstantaneousControl.h
     *  \brief 
     *  
     *  This class derives from EquationSystem and expands its functionalities
     *  by re-defining the solve method for instantaneous control systems.
     *  TODO: maybe TimeDependentEquationSystem from another branch is more suitable
     */
    class InstantaneousControl : public EquationSystem
    {
        // ---------------------------------------------------------------------------------------------//  

        public:
            dolfin::Parameters parameters;

            /******************* TYPEDEFS **********************/
            typedef std::tuple <std::string, int, int> PreviousSolutionLinkValue;
            typedef std::pair <LinkKey, PreviousSolutionLinkValue> PreviousSolutionLink;
            

            /******************* CONSTRUCTORS ******************/
            //! Default constructor 
            InstantaneousControl () = delete;
            
            //InstantaneousControl (std::shared_ptr<dolfin::Form> objectiveFunctional, std::shared_ptr<dolfin::GenericFunction> control);
            InstantaneousControl (std::shared_ptr<dolfin::Form> objectiveFunctional, std::shared_ptr<dolfin::Constant> control);
            
            //InstantaneousControl (std::shared_ptr<dolfin::Form> objectiveFunctional, std::shared_ptr<dolfin::Form> projector, std::shared_ptr<dolfin::GenericFunction> control);
            InstantaneousControl (std::shared_ptr<dolfin::Form> objectiveFunctional, std::shared_ptr<dolfin::Form> projector, std::shared_ptr<dolfin::Constant> control);
            
            /******************* DESTRUCTOR *******************/
            //! Default destructor
            virtual ~InstantaneousControl () = default;
            

            /******************** METHODS *********************/
            virtual void addProblem (const std::string& problemName, dcp::AbstractProblem& problem, bool isTimeDependent=true);

            //! Adds link between problems' coefficient and solution at a previous time step [1]
            /*!
             *  This function will add to the stored map \c linksToPreviousSolutions_ a \c std::pair created on the
             *  input arguments and perform the actual linking calling \c linkProblemToPreviousSolution.
             *  \param linkFrom identifies the problem whose parameter (passed as second argument to the function) 
             *  should be linked with the solution of the problem identified by the fourth parameter (\c linkTo)
             *  \param linkedCoefficientName identifies the coefficient to be linked with said solution
             *  \param linkedCoefficientType identifies the type of the coefficient, and will be passed to the function
             *  \c setCoefficient (see \c dcp::AbstractProblem documentation for more details)
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
             *  \c setCoefficient (see \c dcp::AbstractProblem documentation for more details)
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

            //! Check if system time loop is finished. It basically calls the function \c isFinished() on every problem
            //! stored in \c storedProblems_ and checks if the number of problems whose time loop has ended is equal 
            //! to the size of \c storedProblems_
            /*!
             *  \return true if all problems' time loops have ended, false otherwise
             */
            virtual bool isFinished ();
            
            //! Solve all the problems in the order specified by the private member \c solveOrder_. 
            //! The single problems will be solved calling the \c solve method with \c type argument equal to \c "step"
            /*!
             *  \param forceRelinking a boolean flag which, if set to \c true, overrides the current value of protected 
             *  member variable needsLinksScanning_. Default value is \c false
             */
            virtual void solve (const bool& forceRelinking = false) override;
            
            //! Solve the problem corresponding to the name given [1].
            //! The problem will be solved calling the \c solve method with \c type argument equal to \c "step"
            /*!
             *  \param problemName a string identifying the problem to be solved. If no problem with that name
             *  is found, a warning is issued
             *  \param forceRelinking a boolean flag which, if set to \c true, overrides the current value of protected 
             *  member variable needsLinksScanning_. Default value is \c false
             */
            virtual void solve (const std::string& problemName, const bool& forceRelinking = false) override;
            
            //! Solve the problem corresponding to the name given [2]
            //! The problem will be solved calling the \c solve method with \c type argument equal to \c "step"
            /*!
             *  This method is provided only to allow calls like
             *  \code
             *  solve ("foo_problem");
             *  \endcode
             *  In this case, the compiler would in fact otherwise call \c solve \c (const \c bool&) which is the 
             *  best-matching implicit conversion for a parameter of type \c const \c char*. Using this method,
             *  the version of \c solve that takes a \c std::string is called as expected.
             */
            virtual void solve (const char* problemName, const bool& forceRelinking = false) override;

    //! Set mesh-and-dofs manager
    virtual void setMeshManager (const MeshManager<dolfin::ALE> & meshManager);

    //! Get mesh-and-dofs manager
    virtual dcp::MeshManager<> meshManager () const;

        // ---------------------------------------------------------------------------------------------//  

        protected:

            virtual void solve (const std::string& problemName, const bool& forceRelinking, const bool advanceTimeDependentProblem);

            virtual double chooseIncrementStep (const double alpha, const double alphaMin, dolfin::Array<double> & controlValues, const dolfin::Array<double> & controlPoint, const std::string & controlName, const double & projectedAdjoint, const bool& forceRelinking);

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
             *  Thhis map associates a <tt> std::tuple<std::string, std::string, std::string></tt> 
             *  and a <tt>std::pair <std::string, int></tt>, where:
             *  \li the first \c string contains the name of the problem whose coefficient should be linked against 
             *  some other problem's solution
             *  \li the second \c string contains the type of such coefficient, in a form that can be passed to
             *  \c dcp::AbstractProblem::setCoefficients
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

            //! The stored problems that are unsteady
            std::map <std::string, std::shared_ptr <dcp::AbstractProblem> > timeDependentStoredProblems_;

            const MeshManager<dolfin::ALE> * meshManager_;

            std::shared_ptr<dolfin::Form> objectiveFunctional_;

            std::vector<double> objFunValues_;

            std::shared_ptr<dolfin::Form> projector_;

            //std::shared_ptr<dolfin::GenericFunction> control_;
            std::shared_ptr<dolfin::Constant> control_;

            RescaleThroughTime rescale_;
    };
}

#endif
