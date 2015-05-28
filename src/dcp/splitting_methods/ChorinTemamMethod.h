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

#ifndef SRC_SPLITTING_METHODS_CHORINTEMAMMETHOD_H_INCLUDE_GUARD
#define SRC_SPLITTING_METHODS_CHORINTEMAMMETHOD_H_INCLUDE_GUARD

#include <dcp/splitting_methods/NavierStokesSplittingMethod.h>


namespace dcp
{
    /*! \class ChorinTemamMethod ChorinTemamMethod.h
     *  \brief Class that implements the Chorin-Temam splitting method for
     *  time dependent Navier-Stokes equations.
     *
     *  The three problems in which the Navier-Stokes equations
     *  are split are called \c "prediction_problem", \c "correction_problem"
     *  and \c "projection_problem". These  problems must be passed to the
     *  class through the template arguments. 
     *  The templates arguments are:
     *      - \c T_PredictionBilinearForm the bilinear form of the prediction
     *        problem
     *      - \c T_PredictionLinearForm the linear form of the prediction
     *        problem
     *      - \c T_CorrectionBilinearForm the bilinear form of the correction
     *        problem
     *      - \c T_CorrectionLinearForm the linear form of the correction
     *        problem
     *      - \c T_ProjectionBilinearForm the bilinear form of the projection
     *        problem
     *      - \c T_ProjectionLinearForm the linear form of the projection
     *        problem
     *
     *  The Chorin-Temam problem is supposed to be in the following form:
     *  \f[
     *       \left\{
     *       \begin{array}{l}
     *           u^* - \Delta t \, \nu \Delta u^* + \Delta t \, \left( u \cdot \nabla \right) u^*
     *              = u^n + \Delta t \, f^{n+1} \hfill \mbox{ with b.c.} \\
     *           \Delta t \, \Delta p^{n+1} = \nabla \cdot u^* \hfill \mbox{with b.c.} \\
     *           u^{n+1} = u^* - \Delta t \, \nabla p^{n+1}
     *       \end{array}
     *       \right.
     *  \f]
     */
    template <class T_PredictionBilinearForm_,
              class T_PredictionLinearForm_,
              class T_CorrectionBilinearForm_,
              class T_CorrectionLinearForm_,
              class T_ProjectionBilinearForm_,
              class T_ProjectionLinearForm_>
        class ChorinTemamMethod : public NavierStokesSplittingMethod
        {
            // ---------------------------------------------------------------------------------------------//

            public:
                /************************* TYPEDEFS ************************/
                typedef T_PredictionBilinearForm_ T_PredictionBilinearForm;
                typedef T_PredictionLinearForm_   T_PredictionLinearForm;
                typedef T_CorrectionBilinearForm_ T_CorrectionBilinearForm;
                typedef T_CorrectionLinearForm_   T_CorrectionLinearForm;
                typedef T_ProjectionBilinearForm_ T_ProjectionBilinearForm;
                typedef T_ProjectionLinearForm_   T_ProjectionLinearForm;

                /************************* CONSTRUCTORS ********************/
                //!  Constructor [1]
                /*!
                 *  \param velocityFunctionSpace the velocity function space
                 *  \param pressureFunctionSpace the pressure function space
                 *  \param startTime the initial time for the simulation
                 *  \param dt the time step
                 *  \param endTime the final time for the simulation
                 *  \param nu the kinematic viscosity
                 *  \param previousVelocityName the name of the coefficient representing the velocity at the previous
                 *  time step in the ufl file describing the predicion problem. Default: \c "u_old"
                 *  \param intermediateVelocityName the name of the coefficient representing the intermediate velocity
                 *  computed in the first step of the Chorin-Temam method in the ufl files describing the correction and
                 *  projection problems. Default: \c "u_star"
                 *  \param pressureName the name of the coefficient representing the pressure in the projection problem
                 *  in the ufl files describing the prediction and correction problems. Default: \c "p"
                 *  \param dtName the name of the coefficient representing the time step in the ufl files describing the
                 *  problems. Default: \c "dt"
                 *  
                 *  Note that this constructor will create a \c dcp::Time pointer unique to this object
                 */
                ChorinTemamMethod (const dolfin::FunctionSpace& velocityFunctionSpace,
                                   const dolfin::FunctionSpace& pressureFunctionSpace,
                                   const double& startTime,
                                   const double& dt,
                                   const double& endTime,
                                   const dolfin::GenericFunction& nu,
                                   const std::string& previousVelocityName = "u_old",
                                   const std::string& intermediateVelocityName = "u_star",
                                   const std::string& pressureName = "p",
                                   const std::string& dtName = "dt");

                //!  Constructor [2]
                /*!
                 *  \param velocityFunctionSpace the velocity function space
                 *  \param pressureFunctionSpace the pressure function space
                 *  \param time the shared time object to be used
                 *  \param startTime the initial time for the simulation
                 *  \param dt the time step
                 *  \param endTime the final time for the simulation
                 *  \param nu the kinematic viscosity
                 *  \param previousVelocityName the name of the coefficient representing the velocity at the previous
                 *  time step in the ufl file describing the predicion problem. Default: \c "u_old"
                 *  \param intermediateVelocityName the name of the coefficient representing the intermediate velocity
                 *  computed in the first step of the Chorin-Temam method in the ufl files describing the correction and
                 *  projection problems. Default: \c "u_star"
                 *  \param pressureName the name of the coefficient representing the pressure in the projection problem
                 *  in the ufl files describing the prediction and correction problems. Default: \c "p"
                 *  \param dtName the name of the coefficient representing the time step in the ufl files describing the
                 *  problems. Default: \c "dt"
                 */
                ChorinTemamMethod (const dolfin::FunctionSpace& velocityFunctionSpace,
                                   const dolfin::FunctionSpace& pressureFunctionSpace,
                                   const std::shared_ptr<dcp::Time> time,
                                   const double& startTime,
                                   const double& dt,
                                   const double& endTime,
                                   const dolfin::GenericFunction& nu,
                                   const std::string& previousVelocityName = "u_old",
                                   const std::string& intermediateVelocityName = "u_star",
                                   const std::string& pressureName = "p",
                                   const std::string& dtName = "dt");

                //!  Constructor with time dependent external force
                /*!
                 *  Builds the object when a time dependent external force is used. For a steady external force (which
                 *  is basically a coefficient) use the previous constructor and just call the method 
                 *  \c setCoefficient() on the system contained therein (which can be obtained with the method 
                 *  \c system() ).
                 *  
                 *  \param velocityFunctionSpace the velocity function space
                 *  \param pressureFunctionSpace the pressure function space
                 *  \param startTime the initial time for the simulation
                 *  \param dt the time step
                 *  \param endTime the final time for the simulation
                 *  \param nu the kinematic viscosity
                 *  \param externalForce the external load applied to the prediction problem
                 *  \param previousVelocityName the name of the coefficient representing the velocity at the previous
                 *  time step in the ufl file describing the predicion problem. Default: \c "u_old"
                 *  \param intermediateVelocityName the name of the coefficient representing the intermediate velocity
                 *  computed in the first step of the Chorin-Temam method in the ufl files describing the correction and
                 *  projection problems. Default: \c "u_star"
                 *  \param pressureName the name of the coefficient representing the pressure in the projection problem
                 *  in the ufl files describing the prediction and correction problems. Default: \c "p"
                 *  \param externalForceName the name of the coefficient representing the external force in the 
                 *  prediction problem. Default: "f"
                 *  \param dtName the name of the coefficient representing the time step in the ufl files describing the
                 *  problems. Default: \c "dt"
                 *  
                 *  Note that this constructor will create a \c dcp::Time pointer unique to this object
                 */
                ChorinTemamMethod (const dolfin::FunctionSpace& velocityFunctionSpace,
                                   const dolfin::FunctionSpace& pressureFunctionSpace,
                                   const double& startTime,
                                   const double& dt,
                                   const double& endTime,
                                   const dolfin::GenericFunction& nu,
                                   std::shared_ptr<dcp::TimeDependentExpression> externalForce,
                                   const std::string& previousVelocityName = "u_old",
                                   const std::string& intermediateVelocityName = "u_star",
                                   const std::string& pressureName = "p",
                                   const std::string& externalForceName = "f",
                                   const std::string& dtName = "dt");

                //!  Constructor with time dependent external force [2]
                /*!
                 *  Builds the object when a time dependent external force is used. For a steady external force (which
                 *  is basically a coefficient) use the previous constructor and just call the method 
                 *  \c setCoefficient() on the system contained therein (which can be obtained with the method 
                 *  \c system() ).
                 *  
                 *  \param velocityFunctionSpace the velocity function space
                 *  \param pressureFunctionSpace the pressure function space
                 *  \param time the shared time object to be used
                 *  \param startTime the initial time for the simulation
                 *  \param dt the time step
                 *  \param endTime the final time for the simulation
                 *  \param nu the kinematic viscosity
                 *  \param externalForce the external load applied to the prediction problem
                 *  \param previousVelocityName the name of the coefficient representing the velocity at the previous
                 *  time step in the ufl file describing the predicion problem. Default: \c "u_old"
                 *  \param intermediateVelocityName the name of the coefficient representing the intermediate velocity
                 *  computed in the first step of the Chorin-Temam method in the ufl files describing the correction and
                 *  projection problems. Default: \c "u_star"
                 *  \param pressureName the name of the coefficient representing the pressure in the projection problem
                 *  in the ufl files describing the prediction and correction problems. Default: \c "p"
                 *  \param externalForceName the name of the coefficient representing the external force in the 
                 *  prediction problem. Default: "f"
                 *  \param dtName the name of the coefficient representing the time step in the ufl files describing the
                 *  problems. Default: \c "dt"
                 */
                ChorinTemamMethod (const dolfin::FunctionSpace& velocityFunctionSpace,
                                   const dolfin::FunctionSpace& pressureFunctionSpace,
                                   const std::shared_ptr<dcp::Time> time,
                                   const double& startTime,
                                   const double& dt,
                                   const double& endTime,
                                   const dolfin::GenericFunction& nu,
                                   std::shared_ptr<dcp::TimeDependentExpression> externalForce,
                                   const std::string& previousVelocityName = "u_old",
                                   const std::string& intermediateVelocityName = "u_star",
                                   const std::string& pressureName = "p",
                                   const std::string& externalForceName = "f",
                                   const std::string& dtName = "dt");

                /************************* DESTRUCTOR ********************/
                //! Destructor
                /*!
                 *  Default destructor, since members of the class are trivially destructible.
                 */
                virtual ~ChorinTemamMethod () {};

                // ---------------------------------------------------------------------------------------------//

            protected:

                // ---------------------------------------------------------------------------------------------//

            private:

        };



    // ============================================================================================== //
    // ==================================== IMPLEMENTATION ========================================== //
    // ============================================================================================== //

    /************************* CONSTRUCTORS ********************/
    template <class T_PredictionBilinearForm,
              class T_PredictionLinearForm,
              class T_CorrectionBilinearForm,
              class T_CorrectionLinearForm,
              class T_ProjectionBilinearForm,
              class T_ProjectionLinearForm>
        ChorinTemamMethod<T_PredictionBilinearForm,
                          T_PredictionLinearForm,
                          T_CorrectionBilinearForm,
                          T_CorrectionLinearForm,
                          T_ProjectionBilinearForm,
                          T_ProjectionLinearForm>::
        ChorinTemamMethod (const dolfin::FunctionSpace& velocityFunctionSpace,
                           const dolfin::FunctionSpace& pressureFunctionSpace,
                           const double& startTime,
                           const double& dt,
                           const double& endTime,
                           const dolfin::GenericFunction& nu,
                           const std::string& previousVelocityName,
                           const std::string& intermediateVelocityName,
                           const std::string& pressureName,
                           const std::string& dtName) :
                NavierStokesSplittingMethod 
                    (std::vector<std::shared_ptr <dolfin::FunctionSpace>> 
                         {std::shared_ptr <dolfin::FunctionSpace> (new dolfin::FunctionSpace (velocityFunctionSpace)), 
                          std::shared_ptr <dolfin::FunctionSpace> (new dolfin::FunctionSpace (pressureFunctionSpace))}
                    )
    {
        dolfin::begin (dolfin::DBG, "Building ChorinTemamMethod...");

        parameters ["splitting_method_type"] = "chorin_temam";

        dolfin::log (dolfin::DBG, "Creating the time stepping linear problems...");

        // define chorin temam system time
        std::shared_ptr<dcp::Time> time (new dcp::Time (startTime));

        // define the problems
        // 1) prediction problem
        dolfin::begin (dolfin::DBG, "Creating prediction problem...");

        std::shared_ptr <dcp::AbstractProblem> timeSteppingPredictionProblem
            (new dcp::LinearProblem <T_PredictionBilinearForm, T_PredictionLinearForm> (velocityFunctionSpace_));

        std::shared_ptr <dcp::TimeDependentProblem> predictionProblem
            (new dcp::TimeDependentProblem (timeSteppingPredictionProblem,
                                            time,
                                            startTime,
                                            dt,
                                            endTime,
                                            {"bilinear_form"},
                                            {}));

        predictionProblem->parameters ["dt_name"] = dtName;
        predictionProblem->parameters ["previous_solution_name"] = previousVelocityName;
        predictionProblem->parameters ["previous_solution_is_set_externally"] = true;

        dolfin::end (); // "Creating prediction problem..."


        // 2) correction problem
        dolfin::begin (dolfin::DBG, "Creating correction problem...");

        std::shared_ptr <dcp::AbstractProblem> timeSteppingCorrectionProblem
            (new dcp::LinearProblem <T_CorrectionBilinearForm, T_CorrectionLinearForm> (pressureFunctionSpace_));

        std::shared_ptr <dcp::TimeDependentProblem> correctionProblem
            (new dcp::TimeDependentProblem (timeSteppingCorrectionProblem,
                                            time,
                                            startTime,
                                            dt,
                                            endTime,
                                            {"bilinear_form"},
                                            {}));

        correctionProblem->parameters ["dt_name"] = dtName;

        dolfin::end (); // "Creating correction problem..."


        // 3) projection problem
        dolfin::begin (dolfin::DBG, "Creating projection problem...");

        std::shared_ptr <dcp::AbstractProblem> timeSteppingProjectionProblem
            (new dcp::LinearProblem <T_ProjectionBilinearForm, T_ProjectionLinearForm> (velocityFunctionSpace_));

        std::shared_ptr <dcp::TimeDependentProblem> projectionProblem
            (new dcp::TimeDependentProblem (timeSteppingProjectionProblem,
                                            time,
                                            startTime,
                                            dt,
                                            endTime,
                                            {"linear_form"},
                                            {}));

        projectionProblem->parameters ["dt_name"] = dtName;


        dolfin::end (); // "Creating projection problem..."


        // define the system
        dolfin::begin ("Creating time dependent Chorin-Temam system...");

        // 0) create the object
        std::shared_ptr<dcp::TimeDependentEquationSystem> chorinTemamSystem 
            (new dcp::TimeDependentEquationSystem (time, startTime, dt, endTime));

        // 1) add problems
        dolfin::begin (dolfin::DBG, "Adding problems to protected member map...");
        chorinTemamSystem->addProblem ("prediction_problem", predictionProblem);
        chorinTemamSystem->addProblem ("correction_problem", correctionProblem);
        chorinTemamSystem->addProblem ("projection_problem", projectionProblem);
        dolfin::end (); // "Adding problems to protected member map"

        // 2) add links
        dolfin::begin (dolfin::DBG, "Setting up problems' links...");
        chorinTemamSystem->addLink ("prediction_problem", previousVelocityName, "bilinear_form", "projection_problem");
        chorinTemamSystem->addLink ("prediction_problem", previousVelocityName, "linear_form", "projection_problem");
        chorinTemamSystem->addLink ("correction_problem", intermediateVelocityName, "linear_form", "prediction_problem");
        chorinTemamSystem->addLink ("projection_problem", intermediateVelocityName, "linear_form", "prediction_problem");
        chorinTemamSystem->addLink ("projection_problem", pressureName, "linear_form", "correction_problem");
        dolfin::end (); // "Setting up problems' links"

        // 3) set coefficients
        dolfin::begin (dolfin::DBG, "Setting coefficients...");
        (*chorinTemamSystem) ["prediction_problem"].setCoefficient ("bilinear_form",
                                                                    dolfin::reference_to_no_delete_pointer (nu),
                                                                    "nu");
        dolfin::end (); // "Setting coefficients"

        dolfin::end (); // "Creating the time stepping linear problems..."

        dolfin::begin (dolfin::DBG, "Saving time dependent problem as protected member...");
        differentialSystem_ = chorinTemamSystem;
        dolfin::end (); // "Saving time dependent problem as protected member..."

        dolfin::end (); // "Building ChorinTemamMethod"

        dolfin::log (dolfin::DBG, "ChorinTemamMethod object created");
    }
    


    template <class T_PredictionBilinearForm,
              class T_PredictionLinearForm,
              class T_CorrectionBilinearForm,
              class T_CorrectionLinearForm,
              class T_ProjectionBilinearForm,
              class T_ProjectionLinearForm>
        ChorinTemamMethod<T_PredictionBilinearForm,
                          T_PredictionLinearForm,
                          T_CorrectionBilinearForm,
                          T_CorrectionLinearForm,
                          T_ProjectionBilinearForm,
                          T_ProjectionLinearForm>::
        ChorinTemamMethod (const dolfin::FunctionSpace& velocityFunctionSpace,
                           const dolfin::FunctionSpace& pressureFunctionSpace,
                           const std::shared_ptr<dcp::Time> time,
                           const double& startTime,
                           const double& dt,
                           const double& endTime,
                           const dolfin::GenericFunction& nu,
                           const std::string& previousVelocityName,
                           const std::string& intermediateVelocityName,
                           const std::string& pressureName,
                           const std::string& dtName) :
                NavierStokesSplittingMethod 
                    (std::vector<std::shared_ptr <dolfin::FunctionSpace>> 
                         {std::shared_ptr <dolfin::FunctionSpace> (new dolfin::FunctionSpace (velocityFunctionSpace)), 
                          std::shared_ptr <dolfin::FunctionSpace> (new dolfin::FunctionSpace (pressureFunctionSpace))}
                    )
    {
        dolfin::begin (dolfin::DBG, "Building ChorinTemamMethod...");

        parameters ["splitting_method_type"] = "chorin_temam";

        dolfin::log (dolfin::DBG, "Creating the time stepping linear problems...");


        // define the problems
        // 1) prediction problem
        dolfin::begin (dolfin::DBG, "Creating prediction problem...");

        std::shared_ptr <dcp::AbstractProblem> timeSteppingPredictionProblem
            (new dcp::LinearProblem <T_PredictionBilinearForm, T_PredictionLinearForm> (velocityFunctionSpace_));

        std::shared_ptr <dcp::TimeDependentProblem> predictionProblem
            (new dcp::TimeDependentProblem (timeSteppingPredictionProblem,
                                            time,
                                            startTime,
                                            dt,
                                            endTime,
                                            {"bilinear_form"},
                                            {}));

        predictionProblem->parameters ["dt_name"] = dtName;
        predictionProblem->parameters ["previous_solution_name"] = previousVelocityName;
        predictionProblem->parameters ["previous_solution_is_set_externally"] = true;

        dolfin::end (); // "Creating prediction problem..."


        // 2) correction problem
        dolfin::begin (dolfin::DBG, "Creating correction problem...");

        std::shared_ptr <dcp::AbstractProblem> timeSteppingCorrectionProblem
            (new dcp::LinearProblem <T_CorrectionBilinearForm, T_CorrectionLinearForm> (pressureFunctionSpace_));

        std::shared_ptr <dcp::TimeDependentProblem> correctionProblem
            (new dcp::TimeDependentProblem (timeSteppingCorrectionProblem,
                                            time,
                                            startTime,
                                            dt,
                                            endTime,
                                            {"bilinear_form"},
                                            {}));

        correctionProblem->parameters ["dt_name"] = dtName;

        dolfin::end (); // "Creating correction problem..."


        // 3) projection problem
        dolfin::begin (dolfin::DBG, "Creating projection problem...");

        std::shared_ptr <dcp::AbstractProblem> timeSteppingProjectionProblem
            (new dcp::LinearProblem <T_ProjectionBilinearForm, T_ProjectionLinearForm> (velocityFunctionSpace_));

        std::shared_ptr <dcp::TimeDependentProblem> projectionProblem
            (new dcp::TimeDependentProblem (timeSteppingProjectionProblem,
                                            time, 
                                            startTime,
                                            dt,
                                            endTime,
                                            {"linear_form"},
                                            {}));

        projectionProblem->parameters ["dt_name"] = dtName;


        dolfin::end (); // "Creating projection problem..."


        // define the system
        dolfin::begin ("Creating time dependent Chorin-Temam system...");

        // 0) create the object
        std::shared_ptr<dcp::TimeDependentEquationSystem> chorinTemamSystem 
            (new dcp::TimeDependentEquationSystem (time, startTime, dt, endTime));

        // 1) add problems
        dolfin::begin (dolfin::DBG, "Adding problems to protected member map...");
        chorinTemamSystem->addProblem ("prediction_problem", predictionProblem);
        chorinTemamSystem->addProblem ("correction_problem", correctionProblem);
        chorinTemamSystem->addProblem ("projection_problem", projectionProblem);
        dolfin::end (); // "Adding problems to protected member map"

        // 2) add links
        dolfin::begin (dolfin::DBG, "Setting up problems' links...");
        chorinTemamSystem->addLink ("prediction_problem", previousVelocityName, "bilinear_form", "projection_problem");
        chorinTemamSystem->addLink ("prediction_problem", previousVelocityName, "linear_form", "projection_problem");
        chorinTemamSystem->addLink ("correction_problem", intermediateVelocityName, "linear_form", "prediction_problem");
        chorinTemamSystem->addLink ("projection_problem", intermediateVelocityName, "linear_form", "prediction_problem");
        chorinTemamSystem->addLink ("projection_problem", pressureName, "linear_form", "correction_problem");
        dolfin::end (); // "Setting up problems' links"

        // 3) set coefficients
        dolfin::begin (dolfin::DBG, "Setting coefficients...");
        (*chorinTemamSystem) ["prediction_problem"].setCoefficient ("bilinear_form",
                                                                    dolfin::reference_to_no_delete_pointer (nu),
                                                                    "nu");
        dolfin::end (); // "Setting coefficients"

        dolfin::end (); // "Creating the time stepping linear problems..."

        dolfin::begin (dolfin::DBG, "Saving time dependent problem as protected member...");
        differentialSystem_ = chorinTemamSystem;
        dolfin::end (); // "Saving time dependent problem as protected member..."

        dolfin::end (); // "Building ChorinTemamMethod"

        dolfin::log (dolfin::DBG, "ChorinTemamMethod object created");
    }
    


    template <class T_PredictionBilinearForm,
              class T_PredictionLinearForm,
              class T_CorrectionBilinearForm,
              class T_CorrectionLinearForm,
              class T_ProjectionBilinearForm,
              class T_ProjectionLinearForm>
        ChorinTemamMethod<T_PredictionBilinearForm,
                          T_PredictionLinearForm,
                          T_CorrectionBilinearForm,
                          T_CorrectionLinearForm,
                          T_ProjectionBilinearForm,
                          T_ProjectionLinearForm>::
        ChorinTemamMethod (const dolfin::FunctionSpace& velocityFunctionSpace,
                           const dolfin::FunctionSpace& pressureFunctionSpace,                      
                           const double& startTime,
                           const double& dt,
                           const double& endTime,
                           const dolfin::GenericFunction& nu,
                           std::shared_ptr<dcp::TimeDependentExpression> externalForce,
                           const std::string& previousVelocityName,
                           const std::string& intermediateVelocityName,
                           const std::string& pressureName,
                           const std::string& externalForceName,
                           const std::string& dtName) :
                NavierStokesSplittingMethod 
                    (std::vector<std::shared_ptr <dolfin::FunctionSpace>> 
                         {std::shared_ptr <dolfin::FunctionSpace> (new dolfin::FunctionSpace (velocityFunctionSpace)), 
                          std::shared_ptr <dolfin::FunctionSpace> (new dolfin::FunctionSpace (pressureFunctionSpace))}
                    )
    {
        ChorinTemamMethod (velocityFunctionSpace,
                           pressureFunctionSpace,
                           startTime,
                           dt,
                           endTime,
                           nu,
                           previousVelocityName,
                           intermediateVelocityName,
                           pressureName,
                           dtName);
                          
        (*differentialSystem_) ["prediction_problem"].addTimeDependentCoefficient (externalForceName, "linear_form", externalForce);
    }
    


    template <class T_PredictionBilinearForm,
              class T_PredictionLinearForm,
              class T_CorrectionBilinearForm,
              class T_CorrectionLinearForm,
              class T_ProjectionBilinearForm,
              class T_ProjectionLinearForm>
        ChorinTemamMethod<T_PredictionBilinearForm,
                          T_PredictionLinearForm,
                          T_CorrectionBilinearForm,
                          T_CorrectionLinearForm,
                          T_ProjectionBilinearForm,
                          T_ProjectionLinearForm>::
        ChorinTemamMethod (const dolfin::FunctionSpace& velocityFunctionSpace,
                           const dolfin::FunctionSpace& pressureFunctionSpace,                      
                           const std::shared_ptr<dcp::Time> time,
                           const double& startTime,
                           const double& dt,
                           const double& endTime,
                           const dolfin::GenericFunction& nu,
                           std::shared_ptr<dcp::TimeDependentExpression> externalForce,
                           const std::string& previousVelocityName,
                           const std::string& intermediateVelocityName,
                           const std::string& pressureName,
                           const std::string& externalForceName,
                           const std::string& dtName) :
                NavierStokesSplittingMethod 
                    (std::vector<std::shared_ptr <dolfin::FunctionSpace>> 
                         {std::shared_ptr <dolfin::FunctionSpace> (new dolfin::FunctionSpace (velocityFunctionSpace)), 
                          std::shared_ptr <dolfin::FunctionSpace> (new dolfin::FunctionSpace (pressureFunctionSpace))}
                    )
    {
        ChorinTemamMethod (velocityFunctionSpace,
                           pressureFunctionSpace,
                           time,
                           startTime,
                           dt,
                           endTime,
                           nu,
                           previousVelocityName,
                           intermediateVelocityName,
                           pressureName,
                           dtName);
                          
        (*differentialSystem_) ["prediction_problem"].addTimeDependentCoefficient (externalForceName, "linear_form", externalForce);
    }
}
#endif

