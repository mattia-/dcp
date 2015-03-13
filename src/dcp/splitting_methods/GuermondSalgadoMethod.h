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

#ifndef SRC_SPLITTING_METHODS_GUERMONDSALGADOMETHOD_H_INCLUDE_GUARD
#define SRC_SPLITTING_METHODS_GUERMONDSALGADOMETHOD_H_INCLUDE_GUARD

#include <dcp/splitting_methods/NavierStokesSplittingMethod.h>


namespace dcp
{
    /*! \class GuermondSalgadoMethod GuermondSalgadoMethod.h
     *  \brief Class that implements the incremental E Guermond-Salgado splitting
     *  method for time dependent Navier-Stokes equations with variable density,
     *  as described in 
     *  "J.-L. Guermond, A. Salgado, A splitting method for incompressible flows with variable
     *  density ..., J. Comput. Phys. (2009), doi:10.1016/j.jcp.2008.12.036."
     *
     *  The four problems in which the Navier-Stokes equations
     *  are split are called \c "density_problem", \c "velocity_problem",
     *  \c "pressure_correction_problem" and \c "pressure_update_problem" (with 
     *  self-explanatory names).
     *  These problems must be passed to the class through the template arguments.
     *  The templates arguments are:
     *      - \c T_DensityBilinearForm the bilinear form of the density
     *        problem
     *      - \c T_DensityLinearForm the linear form of the density
     *        problem
     *      - \c T_VelocityBilinearForm the bilinear form of the velocity
     *        problem
     *      - \c T_VelocityLinearForm the linear form of the velocity
     *        problem
     *      - \c T_PressureCorrectionBilinearForm the bilinear form of the pressure
     *        correction problem
     *      - \c T_PressureCorrectionLinearForm the linear form of the pressure
     *        correction problem
     *      - \c T_PressureUpdateBilinearForm the bilinear form of the pressure
     *        update problem
     *      - \c T_PressureUpdateLinearForm the linear form of the pressure
     *        update problem
     *
     *  The incremental Guermond-Salgado problem is supposed to be in the form presented in section
     *  3.1 (for order one in time) or 4 (for order two in time) of the aforementioned paper.
     */
    template <class T_DensityBilinearForm_,
              class T_DensityLinearForm_,
              class T_VelocityBilinearForm_,
              class T_VelocityLinearForm_,
              class T_PressureCorrectionBilinearForm_,
              class T_PressureCorrectionLinearForm_,
              class T_PressureUpdateBilinearForm_,
              class T_PressureUpdateLinearForm_>
        class GuermondSalgadoMethod : public NavierStokesSplittingMethod
        {
            // ---------------------------------------------------------------------------------------------//

            public:
                /************************* TYPEDEFS ************************/
                typedef T_DensityBilinearForm_            T_DensityBilinearForm;
                typedef T_DensityLinearForm_              T_DensityLinearForm;
                typedef T_VelocityBilinearForm_           T_VelocityBilinearForm;
                typedef T_VelocityLinearForm_             T_VelocityLinearForm;
                typedef T_PressureCorrectionBilinearForm_ T_PressureCorrectionBilinearForm;
                typedef T_PressureCorrectionLinearForm_   T_PressureCorrectionLinearForm;
                typedef T_PressureUpdateBilinearForm_     T_PressureUpdateBilinearForm;
                typedef T_PressureUpdateLinearForm_       T_PressureUpdateLinearForm;

                /************************* CONSTRUCTORS ********************/
                //!  Constructor
                /*!
                 *  \param densityFunctionSpace the density function space
                 *  \param velocityFunctionSpace the velocity function space
                 *  \param pressureFunctionSpace the pressure function space
                 *  \param startTime the initial time for the simulation
                 *  \param dt the time step
                 *  \param endTime the final time for the simulation
                 *  \param mu the viscosity
                 *  \param chi the penalization parameter in the third equation
                 *  \param nTimeSchemeSteps the number of time steps involved in the time stepping problem solution. 
                 *  See \c dcp::TimeDependentProblem documentation for more details. Default value: \c 1
                 *  \param velocityName the name of the coefficient representing the velocity at the current
                 *  time step in the ufl file describing the pressure correction problem
                 *  Default: \c "u"
                 *  \param previousVelocityName the name of the coefficient representing the velocity at the previous
                 *  time steps in the ufl file describing the density problem and the velocity problem. 
                 *  A number will be appended to this variable to identify solutions from 2 or more time steps back.
                 *  For example, if the current time index is \c n and \c previousVelocityName contains \c "u_old", 
                 *  then the velocity at index \c n-1 will be identified by \c "u_old", the velocity at time index
                 *  \c n-1 by "u_old_2" and so on.
                 *  Default: \c "u_old"
                 *  \param previousPressureName the name of the coefficient representing the pressure at the previous
                 *  time step in the ufl file describing the velocity problem and the pressure correction problem. 
                 *  A number will be appended to this variable to identify solutions from 2 or more time steps back.
                 *  For example, if the current time index is \c n and \c previousPressureName contains \c "p_old", 
                 *  then the pressure at index \c n-1 will be identified by \c "p_old", the pressure at time index
                 *  \c n-1 by "p_old_2" and so on.
                 *  Default: \c "p_old"
                 *  \param densityName the name of the coefficient representing the density at the current time step
                 *  in the ufl file describing the velocity problem.
                 *  Default: \c "rho"
                 *  \param previousDensityName the name of the coefficient representing the density at the previous
                 *  time step in the ufl file describing the density problem and the velocity problem. 
                 *  A number will be appended to this variable to identify solutions from 2 or more time steps back.
                 *  For example, if the current time index is \c n and \c previousDensityName contains \c "rho_old", 
                 *  then the density at index \c n-1 will be identified by \c "rho_old", the density at time index
                 *  \c n-1 by "rho_old_2" and so on.
                 *  Default: \c "rho_old"
                 *  \param pressureIncrementName the name of the coefficient representing the pressure increment
                 *  in the pressure update problem. Default: \c "phi"
                 *  \param previousPressureIncrementName the name of the coefficient representing the pressure increment
                 *  at the previous time step in the velocity problem. 
                 *  A number will be appended to this variable to identify solutions from 2 or more time steps back.
                 *  For example, if the current time index is \c n and \c previousPressureIncrementName contains 
                 *  \c "phi_old", then the pressure increment at index \c n-1 will be identified by \c "phi_old", 
                 *  the pressure increment at time index \c n-1 by "phi_old_2" and so on.
                 *  Default: \c "phi_old"
                 *  \param dtName the name of the coefficient representing the time step in the ufl files describing the
                 *  problems. Default: \c "dt"   
                 */
                GuermondSalgadoMethod (const dolfin::FunctionSpace& densityFunctionSpace,
                                       const dolfin::FunctionSpace& velocityFunctionSpace,
                                       const dolfin::FunctionSpace& pressureFunctionSpace,
                                       const double& startTime,
                                       const double& dt,
                                       const double& endTime,
                                       const dolfin::GenericFunction& mu,
                                       const dolfin::GenericFunction& chi,
                                       const unsigned int& nTimeSchemeSteps = 1,
                                       const std::string& velocityName = "u",
                                       const std::string& previousVelocityName = "u_old",
                                       const std::string& previousPressureName = "p_old",
                                       const std::string& densityName = "rho",
                                       const std::string& previousDensityName = "rho_old",
                                       const std::string& pressureIncrementName = "phi",
                                       const std::string& previousPressureIncrementName = "phi_old",
                                       const std::string& dtName = "dt");

                //!  Constructor with time dependent external force
                /*!
                 *  Builds the object when a time dependent external force is used. For a steady external force (which
                 *  is basically a coefficient) use the previous constructor and just call the method 
                 *  \c setCoefficient() on the system contained therein (which can be obtained with the method 
                 *  \c system() ).
                 *  
                 *  \param densityFunctionSpace the density function space
                 *  \param velocityFunctionSpace the velocity function space
                 *  \param pressureFunctionSpace the pressure function space
                 *  \param startTime the initial time for the simulation
                 *  \param dt the time step
                 *  \param endTime the final time for the simulation
                 *  \param mu the viscosity
                 *  \param chi the penalization parameter in the third equation
                 *  \param externalForce the external load applied to the velocity problem
                 *  \param nTimeSchemeSteps the number of time steps involved in the time stepping problem solution. 
                 *  See \c dcp::TimeDependentProblem documentation for more details. Default value: \c 1
                 *  \param velocityName the name of the coefficient representing the velocity at the current
                 *  time step in the ufl file describing the pressure correction problem
                 *  Default: \c "u"
                 *  \param previousVelocityName the name of the coefficient representing the velocity at the previous
                 *  time step in the ufl file describing the density problem and the velocity problem. 
                 *  Default: \c "u_old"
                 *  \param previousPressureName the name of the coefficient representing the pressure at the previous
                 *  time step in the ufl file describing the velocity problem and the pressure correction problem. 
                 *  Default: \c "p_old"
                 *  \param densityName the name of the coefficient representing the density at the current time step
                 *  in the ufl file describing the velocity problem.
                 *  Default: \c "rho"
                 *  \param previousDensityName the name of the coefficient representing the density at the previous
                 *  time step in the ufl file describing the density problem and the velocity problem. 
                 *  Default: \c "rho_old"
                 *  \param pressureIncrementName the name of the coefficient representing the pressure increment
                 *  in the pressure update problem. Default: \c "phi"
                 *  \param previousPressureIncrementName the name of the coefficient representing the pressure increment
                 *  at the previous time step in the velocity problem. Default: \c "phi_old"
                 *  \param externalForceName the name of the coefficient representing the external force in the 
                 *  prediction problem. Default: "f"
                 *  \param dtName the name of the coefficient representing the time step in the ufl files describing the
                 *  problems. Default: \c "dt"   
                 */
                GuermondSalgadoMethod (const dolfin::FunctionSpace& densityFunctionSpace,
                                       const dolfin::FunctionSpace& velocityFunctionSpace,
                                       const dolfin::FunctionSpace& pressureFunctionSpace,
                                       const double& startTime,
                                       const double& dt,
                                       const double& endTime,
                                       const dolfin::GenericFunction& mu,
                                       const dolfin::GenericFunction& chi,
                                       std::shared_ptr<dcp::TimeDependentExpression> externalForce,
                                       const unsigned int& nTimeSchemeSteps = 1,
                                       const std::string& velocityName = "u",
                                       const std::string& previousVelocityName = "u_old",
                                       const std::string& previousPressureName = "p_old",
                                       const std::string& densityName = "rho",
                                       const std::string& previousDensityName = "rho_old",
                                       const std::string& pressureIncrementName = "phi",
                                       const std::string& previousPressureIncrementName = "phi_old",
                                       const std::string& externalForceName = "f",
                                       const std::string& dtName = "dt");

                /************************* DESTRUCTOR *************************/
                //! Destructor
                /*!
                 hi*  Default destructor, since members of the class are trivially destructible.
                 */
                virtual ~GuermondSalgadoMethod () {};

                // ---------------------------------------------------------------------------------------------//

            protected:
                //! The density function space. It is just a reference to the element pointed by the third element of 
                //! \c functionSpaces_, which is inherited from \c dcp::AbstractSplittingMethod.
                const dolfin::FunctionSpace& densityFunctionSpace_;
                
                // ---------------------------------------------------------------------------------------------//

            private:

        };



    // ============================================================================================== //
    // ==================================== IMPLEMENTATION ========================================== //
    // ============================================================================================== //

    /************************* CONSTRUCTORS ********************/
    template <class T_DensityBilinearForm,
              class T_DensityLinearForm,
              class T_VelocityBilinearForm,
              class T_VelocityLinearForm,
              class T_PressureCorrectionBilinearForm,
              class T_PressureCorrectionLinearForm,
              class T_PressureUpdateBilinearForm,
              class T_PressureUpdateLinearForm>
        GuermondSalgadoMethod<T_DensityBilinearForm,
                              T_DensityLinearForm,
                              T_VelocityBilinearForm,
                              T_VelocityLinearForm,
                              T_PressureCorrectionBilinearForm,
                              T_PressureCorrectionLinearForm,
                              T_PressureUpdateBilinearForm,
                              T_PressureUpdateLinearForm>::
        GuermondSalgadoMethod (const dolfin::FunctionSpace& densityFunctionSpace,
                               const dolfin::FunctionSpace& velocityFunctionSpace,
                               const dolfin::FunctionSpace& pressureFunctionSpace,
                               const double& startTime,
                               const double& dt,
                               const double& endTime,
                               const dolfin::GenericFunction& mu,
                               const dolfin::GenericFunction& chi,
                               const unsigned int& nTimeSchemeSteps,
                               const std::string& velocityName,
                               const std::string& previousVelocityName,
                               const std::string& previousPressureName,
                               const std::string& densityName,
                               const std::string& previousDensityName,
                               const std::string& pressureIncrementName,
                               const std::string& previousPressureIncrementName,
                               const std::string& dtName) :
                NavierStokesSplittingMethod 
                    (std::vector<std::shared_ptr <dolfin::FunctionSpace>> 
                         {std::shared_ptr<dolfin::FunctionSpace> (new dolfin::FunctionSpace (velocityFunctionSpace)), 
                          std::shared_ptr<dolfin::FunctionSpace> (new dolfin::FunctionSpace (pressureFunctionSpace)),
                          std::shared_ptr<dolfin::FunctionSpace> (new dolfin::FunctionSpace (densityFunctionSpace))}
                    ),
                densityFunctionSpace_ (*(functionSpaces_ [2]))
    {
        dolfin::begin (dolfin::DBG, "Building GuermondSalgadoMethod...");
        
        parameters ["splitting_method_type"] = "guermond_salgado";
        
        if (nTimeSchemeSteps != 1 && nTimeSchemeSteps != 2)
        {
            dolfin::dolfin_error ("GuermondSalgadoMethod.h",
                                  "create GuermondSalgadoMethod object",
                                  "Unknown order %d for time scheme. Available orders are 1 and 2", 
                                  nTimeSchemeSteps);
        }
        
        dolfin::begin (dolfin::DBG, "Creating the time stepping linear problems...");

        // define the problems
        // 1) density problem
        dolfin::begin (dolfin::DBG, "Creating density problem...");

        std::shared_ptr <dcp::AbstractProblem> timeSteppingDensityProblem
            (new dcp::LinearProblem <T_DensityBilinearForm, T_DensityLinearForm> (densityFunctionSpace_));

        std::shared_ptr <dcp::TimeDependentProblem> densityProblem
            (new dcp::TimeDependentProblem (timeSteppingDensityProblem,
                                            startTime,
                                            dt,
                                            endTime,
                                            {"bilinear_form", "linear_form"},
                                            {"linear_form"},
                                            nTimeSchemeSteps));

        densityProblem->parameters ["dt_name"] = dtName;
        densityProblem->parameters ["previous_solution_name"] = previousDensityName;

        dolfin::end (); // "Creating density problem..."


        // 2) velocity problem
        dolfin::begin (dolfin::DBG, "Creating velocity problem...");

        std::shared_ptr <dcp::AbstractProblem> timeSteppingVelocityProblem
            (new dcp::LinearProblem <T_VelocityBilinearForm, T_VelocityLinearForm> (velocityFunctionSpace_));

        std::shared_ptr <dcp::TimeDependentProblem> velocityProblem
            (new dcp::TimeDependentProblem (timeSteppingVelocityProblem,
                                            startTime,
                                            dt,
                                            endTime,
                                            {"bilinear_form", "linear_form"},
                                            {"bilinear_form", "linear_form"},
                                            nTimeSchemeSteps));

        velocityProblem->parameters ["dt_name"] = dtName;
        velocityProblem->parameters ["previous_solution_name"] = previousVelocityName;

        dolfin::end (); // "Creating correction problem..."


        // 3) pressure correction problem
        dolfin::begin (dolfin::DBG, "Creating pressure correction problem...");

        std::shared_ptr <dcp::AbstractProblem> timeSteppingPressureCorrectionProblem
            (new dcp::LinearProblem <T_PressureCorrectionBilinearForm, 
             T_PressureCorrectionLinearForm> 
             (pressureFunctionSpace_));

        std::shared_ptr <dcp::TimeDependentProblem> pressureCorrectionProblem
            (new dcp::TimeDependentProblem (timeSteppingPressureCorrectionProblem,
                                            startTime,
                                            dt,
                                            endTime,
                                            {"linear_form"},
                                            {},
                                            nTimeSchemeSteps));

        pressureCorrectionProblem->parameters ["dt_name"] = dtName;


        dolfin::end (); // "Creating pressure correction problem..."


        // 4) pressure update problem
        dolfin::begin (dolfin::DBG, "Creating pressure update problem...");

        std::shared_ptr <dcp::AbstractProblem> timeSteppingPressureUpdateProblem
            (new dcp::LinearProblem <T_PressureUpdateBilinearForm, 
             T_PressureUpdateLinearForm> 
             (pressureFunctionSpace_));

        // if nTimeSchemeSteps == 2, the problem pressure update problem is created with a start time of startTime+dt, 
        // because otherwise the time would be skewed with respect to the other problems, since they all use two-steps 
        // time schemes (see also dcp::TimeDependentProblem constructor to see why multi-step problems have an 
        // incremented start time)
        std::shared_ptr <dcp::TimeDependentProblem> pressureUpdateProblem = nullptr;
        if (nTimeSchemeSteps == 1)
        {
            pressureUpdateProblem.reset (new dcp::TimeDependentProblem (timeSteppingPressureUpdateProblem,
                                                                        startTime,
                                                                        dt,
                                                                        endTime,
                                                                        {},
                                                                        {"linear_form"}));
        }
        else
        {
            pressureUpdateProblem.reset (new dcp::TimeDependentProblem (timeSteppingPressureUpdateProblem,
                                                                        startTime+dt,
                                                                        dt,
                                                                        endTime,
                                                                        {},
                                                                        {"linear_form"}));
        }

        pressureUpdateProblem->parameters ["previous_solution_name"] = previousPressureName;

        dolfin::end (); // "Creating pressure update problem..."

        dolfin::end (); // "Creating the time stepping linear problems..."

        // define the system
        dolfin::begin (dolfin::DBG, "Creating time dependent Guermond-Salgado system...");

        // 0) create the object
        std::shared_ptr<dcp::TimeDependentEquationSystem> guermondSalgadoSystem (new dcp::TimeDependentEquationSystem);

        // 1) add problems
        dolfin::begin (dolfin::DBG, "Adding problems to protected member map...");
        guermondSalgadoSystem->addProblem ("density_problem", densityProblem);
        guermondSalgadoSystem->addProblem ("velocity_problem", velocityProblem);
        guermondSalgadoSystem->addProblem ("pressure_correction_problem", pressureCorrectionProblem);
        guermondSalgadoSystem->addProblem ("pressure_update_problem", pressureUpdateProblem);
        dolfin::end (); // "Adding problems to protected member map"

        // 2) add links
        dolfin::begin (dolfin::DBG, "Setting up problems' links...");
        guermondSalgadoSystem->addLink ("density_problem", 
                                        previousVelocityName, 
                                        "bilinear_form", 
                                        "velocity_problem");
        guermondSalgadoSystem->addLink ("velocity_problem", 
                                        densityName, 
                                        "bilinear_form", 
                                        "density_problem");
        guermondSalgadoSystem->addLink ("velocity_problem", 
                                        previousPressureName, 
                                        "linear_form", 
                                        "pressure_update_problem");
        guermondSalgadoSystem->addLink ("velocity_problem", 
                                        previousPressureIncrementName, 
                                        "linear_form", 
                                        "pressure_correction_problem");
        guermondSalgadoSystem->addLink ("pressure_correction_problem", 
                                        velocityName, 
                                        "linear_form", 
                                        "velocity_problem");
        guermondSalgadoSystem->addLink ("pressure_update_problem", 
                                        pressureIncrementName, 
                                        "linear_form", 
                                        "pressure_correction_problem");
        
        dolfin::end (); // "Setting up problems' links"
        
        // add links specific to first order time scheme
        if (nTimeSchemeSteps == 1)
        {
            dolfin::begin (dolfin::DBG, "Setting up specific links for first order time scheme...");
            guermondSalgadoSystem->addLinkToPreviousSolution ("velocity_problem", 
                                                              previousDensityName, 
                                                              "bilinear_form", 
                                                              "density_problem",
                                                              1);
            guermondSalgadoSystem->addLinkToPreviousSolution ("velocity_problem", 
                                                              previousDensityName, 
                                                              "linear_form", 
                                                              "density_problem",
                                                              1);
            dolfin::end (); // "Setting up specific links for first order time scheme"
        }
            
        // add links specific to second order time scheme
        if (nTimeSchemeSteps == 2)
        {
            dolfin::begin (dolfin::DBG, "Setting up specific links for second order time scheme...");
            guermondSalgadoSystem->addLinkToPreviousSolution ("density_problem", 
                                                              previousVelocityName + "_2", 
                                                              "bilinear_form", 
                                                              "velocity_problem",
                                                              1);
            guermondSalgadoSystem->addLinkToPreviousSolution ("velocity_problem", 
                                                              previousPressureIncrementName + "_2", 
                                                              "linear_form", 
                                                              "pressure_correction_problem",
                                                              1);
            guermondSalgadoSystem->addLink ("pressure_update_problem", 
                                            velocityName, 
                                            "linear_form", 
                                            "velocity_problem");
            dolfin::end (); // "Setting up specific links for second order time scheme"
        }
            

        // 3) set coefficients
        dolfin::begin (dolfin::DBG, "Setting coefficients...");
        (*guermondSalgadoSystem) ["velocity_problem"].setCoefficient ("bilinear_form",
                                                                      dolfin::reference_to_no_delete_pointer (mu),
                                                                      "mu");
        if (nTimeSchemeSteps == 2)
        {
            (*guermondSalgadoSystem) ["pressure_update_problem"].setCoefficient ("linear_form",
                                                                                 dolfin::reference_to_no_delete_pointer (mu),
                                                                                 "mu");
        }
        (*guermondSalgadoSystem) ["pressure_correction_problem"].setCoefficient 
            ("linear_form",
             dolfin::reference_to_no_delete_pointer (chi),
             "chi");
        dolfin::end (); // "Setting coefficients"

        dolfin::end (); // "Creating the time stepping linear problems..."

        dolfin::begin (dolfin::DBG, "Saving time dependent problem as protected member...");
        differentialSystem_ = guermondSalgadoSystem;
        dolfin::end (); // "Saving time dependent problem as protected member..."

        dolfin::end (); // "Building GuermondSalgadoMethod"

        dolfin::log (dolfin::DBG, "GuermondSalgadoMethod object created");
    }
    
    

    template <class T_DensityBilinearForm,
              class T_DensityLinearForm,
              class T_VelocityBilinearForm,
              class T_VelocityLinearForm,
              class T_PressureCorrectionBilinearForm,
              class T_PressureCorrectionLinearForm,
              class T_PressureUpdateBilinearForm,
              class T_PressureUpdateLinearForm>
        GuermondSalgadoMethod<T_DensityBilinearForm,
                              T_DensityLinearForm,
                              T_VelocityBilinearForm,
                              T_VelocityLinearForm,
                              T_PressureCorrectionBilinearForm,
                              T_PressureCorrectionLinearForm,
                              T_PressureUpdateBilinearForm,
                              T_PressureUpdateLinearForm>::
        GuermondSalgadoMethod (const dolfin::FunctionSpace& densityFunctionSpace,
                               const dolfin::FunctionSpace& velocityFunctionSpace,
                               const dolfin::FunctionSpace& pressureFunctionSpace,
                               const double& startTime,
                               const double& dt,
                               const double& endTime,
                               const dolfin::GenericFunction& mu,
                               const dolfin::GenericFunction& chi,
                               std::shared_ptr<dcp::TimeDependentExpression> externalForce,
                               const unsigned int& nTimeSchemeSteps,
                               const std::string& velocityName,
                               const std::string& previousVelocityName,
                               const std::string& previousPressureName,
                               const std::string& densityName,
                               const std::string& previousDensityName,
                               const std::string& pressureIncrementName,
                               const std::string& previousPressureIncrementName,
                               const std::string& externalForceName,
                               const std::string& dtName) :
                NavierStokesSplittingMethod 
                    (std::vector<std::shared_ptr <dolfin::FunctionSpace>> 
                         {std::shared_ptr<dolfin::FunctionSpace> (new dolfin::FunctionSpace (velocityFunctionSpace)), 
                          std::shared_ptr<dolfin::FunctionSpace> (new dolfin::FunctionSpace (pressureFunctionSpace)),
                          std::shared_ptr<dolfin::FunctionSpace> (new dolfin::FunctionSpace (densityFunctionSpace))}
                    ),
                densityFunctionSpace_ (*(functionSpaces_ [2]))
    {
        GuermondSalgadoMethod (densityFunctionSpace,
                               velocityFunctionSpace,
                               pressureFunctionSpace,
                               startTime,
                               dt,
                               endTime,
                               mu,
                               chi,
                               nTimeSchemeSteps,
                               velocityName,
                               previousVelocityName,
                               previousPressureName,
                               densityName,
                               previousDensityName,
                               pressureIncrementName,
                               previousPressureIncrementName,
                               dtName);
                               
        (*differentialSystem_) ["velocity_problem"].addTimeDependentCoefficient (externalForceName, 
                                                                                 "linear_form", 
                                                                                 externalForce);
    }
}
#endif


