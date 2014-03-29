#ifndef SRC_DIFFERENTIALPROBLEM_NONLINEARDIFFERENTIALPROBLEM_HPP_INCLUDE_GUARD
#define SRC_DIFFERENTIALPROBLEM_NONLINEARDIFFERENTIALPROBLEM_HPP_INCLUDE_GUARD

#include <dolfin.h>
#include <vector>
#include <map>
#include <string>
#include <memory>
#include <DifferentialProblem/AbstractDifferentialProblem.hpp>

namespace control_problem
{
    /*! \class NonlinearDifferentialProblem NonlinearDifferentialProblem.hpp
     *  \brief Class for non-linear differential problems.
     *
     *  This class represents problem of the form
     *  \f[
     *      \mbox{Find } u \in V : F \left(u, v\right) = 0 \ \forall\,v\,\in\,V
     *  \f]
     *  with \f$ F \left(u, v\right) : V \times V \rightarrow \mathds{R}\f$ non linear in the problem unknown \f$u\f$.
     *  
     *  It inherits publicly from \c AbstractDifferentialProblem
     *  and it extends its functionalities to a concrete differential
     *  problem.
     *  Template arguments are:
     *  \arg T_ResidualForm the residual form type, that is the type that describes the form \f$F\f$
     *  \arg T_JacobianForm the jacobian form type, that is the type of the derivative of \f$F\f$ with respect to \f$u\f$,
     *  used in the Newton-Raphson method iterations
     *  This is needed because the dolfin function "coefficient_number" was neither declared
     *  virtual nor implemented by FEniCS's developers in the dolfin class "Form", from which
     *  we first thought to derive for the protected members of the NonlinearDifferentialProblem
     *  class
     */

    template <class T_ResidualForm_, class T_JacobianForm_>
        class NonlinearDifferentialProblem : public AbstractDifferentialProblem
        {
            // ---------------------------------------------------------------------------------------------//  

            public:
                /******************* TYPEDEFS *******************/
                typedef T_ResidualForm_ T_ResidualForm;
                typedef T_JacobianForm_ T_JacobianForm;
                
                
                /******************* CONSTRUCTORS *******************/
                //! Default constructor is deleted. The class is not default constructable.
                NonlinearDifferentialProblem () = delete;

                //!  Constructor with shared pointers [1]
                /*!
                 *  \param mesh the problem mesh as a const std::shared_ptr to dolfin::Mesh
                 *  \param functionSpace the problem finite element space as a const std::shared_ptr to dolfin::FunctionSpace
                 *  \param residualFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the ResidualForm 
                 *  \param jacobianFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the JacobianForm. Default value is the empty string, in which case \c residualFormSolutionName
                 *  will be used for jacobianFormSolutionName
                 *  The stored mesh's and function space's ownership will be shared between the object and the input argument.
                 *  The residual and jacobian form will be created too, calling the constructor which takes the function space
                 *  as input.
                 */
                NonlinearDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
                                              const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
                                              const std::string& residualFormSolutionName,
                                              const std::string& jacobianFormSolutionName);
                

                //! Constructor with references [1]
                /*!
                 *  \param mesh the problem mesh as a const dolfin::Mesh&
                 *  \param functionSpace the problem finite element space as a const dolfin::FunctionSpace&
                 *  \param residualFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the ResidualForm 
                 *  \param jacobianFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the JacobianForm. Default value is the empty string, in which case \c residualFormSolutionName
                 *  will be used for jacobianFormSolutionName
                 *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
                 *  initialized using the \c new operator and mesh's and functionSpace's copy constructor
                 *  The residual and jacobian form will be created too, calling the constructor which takes the function space
                 *  as input.
                 */
                NonlinearDifferentialProblem (const dolfin::Mesh& mesh, 
                                              const dolfin::FunctionSpace& functionSpace,
                                              const std::string& residualFormSolutionName,
                                              const std::string& jacobianFormSolutionName);

                //! Constructor with rvalue references [1]
                /*!
                 *  \param mesh the problem mesh as a dolfin::Mesh&&
                 *  \param functionSpace the problem finite element space as a dolfin::FunctionSpace&&
                 *  \param residualFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the ResidualForm 
                 *  \param jacobianFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the JacobianForm. Default value is the empty string, in which case \c residualFormSolutionName
                 *  will be used for jacobianFormSolutionName
                 *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
                 *  initialized using the \c new operator and mesh's and functionSpace's move constructor
                 *  The residual and jacobian form will be created too, calling the constructor which takes the function space
                 *  as input.
                 */
                NonlinearDifferentialProblem (dolfin::Mesh&& mesh, 
                                              dolfin::FunctionSpace&& functionSpace,
                                              const std::string& residualFormSolutionName,
                                              const std::string& jacobianFormSolutionName);

                
                //!  Constructor with shared pointers [2]
                /*!
                 *  \param mesh the problem mesh as a const std::shared_ptr to dolfin::Mesh
                 *  \param functionSpace the problem finite element space as a const std::shared_ptr to dolfin::FunctionSpace
                 *  \param residualForm a const reference to the problem's residual form
                 *  \param jacobianForm a const reference to the problem's jacobian form
                 *  \param residualFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the ResidualForm 
                 *  \param jacobianFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the JacobianForm. Default value is the empty string, in which case \c residualFormSolutionName
                 *  will be used for jacobianFormSolutionName
                 *  The stored mesh's and function space's ownership will be shared between the object and the input argument.
                 *  The residual and jacobian form will be created too, calling the constructor which takes the function space
                 *  as input.
                 */
                NonlinearDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
                                              const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
                                              const T_ResidualForm& residualForm,
                                              const T_JacobianForm& jacobianForm,
                                              const std::string& residualFormSolutionName,
                                              const std::string& jacobianFormSolutionName);

                //! Constructor with references [2]
                /*!
                 *  \param mesh the problem mesh as a const dolfin::Mesh&
                 *  \param functionSpace the problem finite element space as a const dolfin::FunctionSpace&
                 *  \param residualForm a const reference to the problem's residual form
                 *  \param jacobianForm a const reference to the problem's jacobian form
                 *  \param residualFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the ResidualForm 
                 *  \param jacobianFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the JacobianForm. Default value is the empty string, in which case \c residualFormSolutionName
                 *  will be used
                 *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
                 *  initialized using the \c new operator and mesh's and functionSpace's copy constructor
                 *  The residual and jacobian form will be created too, calling the constructor which takes the function space
                 *  as input.
                 */
                NonlinearDifferentialProblem (const dolfin::Mesh& mesh, 
                                              const dolfin::FunctionSpace& functionSpace,
                                              const T_ResidualForm& residualForm,
                                              const T_JacobianForm& jacobianForm,
                                              const std::string& residualFormSolutionName,
                                              const std::string& jacobianFormSolutionName);

                //! Constructor with rvalue references [3]
                /*!
                 *  \param mesh the problem mesh as a dolfin::Mesh&&
                 *  \param functionSpace the problem finite element space as a dolfin::FunctionSpace&&
                 *  \param residualForm a rvalue reference to the problem's residual form
                 *  \param jacobianForm a rvalue reference to the problem's jacobian form
                 *  \param residualFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the ResidualForm 
                 *  \param jacobianFormSolutionName a string that identifies the name of the function representing the problem
                 *  solution in the JacobianForm. Default value is the empty string, in which case \c residualFormSolutionName
                 *  will be used for jacobianFormSolutionName
                 *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
                 *  initialized using the \c new operator and mesh's and functionSpace's move constructor
                 *  The residual and jacobian form will be created too, calling the constructor which takes the function space
                 *  as input.
                 */
                NonlinearDifferentialProblem (dolfin::Mesh&& mesh, 
                                              dolfin::FunctionSpace&& functionSpace,
                                              T_ResidualForm&& residualForm,
                                              T_JacobianForm&& jacobianForm,
                                              const std::string& residualFormSolutionName,
                                              const std::string& jacobianFormSolutionName);

                /******************* DESTRUCTOR *******************/
                
                //! Destructor
                /*! Default destructor, since members of the class are trivially 
                 * destructible.
                 * It is declared virtual so that derived classes' constructor
                 * can be called on derived classes.
                 * The "default-ness" is set in implementation outside of the class for compatibility with
                 * gcc-4.6, which does not allow virtual members to be defaulted in class
                 */
                virtual ~NonlinearDifferentialProblem ();

                
                /******************* GETTERS *******************/
                //! Get const reference to the problem's residual form
                /*! 
                 *  \return a const reference to the problem's residual form
                 */
                const T_ResidualForm& residualForm () const;

                //! Get const reference to the problem's jacobian form
                /*! 
                 *  \return a const reference to the problem's jacobian form
                 */
                const T_JacobianForm& jacobianForm () const;

                /******************* SETTERS *******************/

                //! Set coefficient [1]. Override of virtual function in \c AbstractDifferentialProblem.
                /*!
                 *  Possible values for \c coefficientType are:
                 *  \li residual_form to set the coefficient in the residual form
                 *  \li jacobian_form to set the coefficient in the jacobian form
                 *  \li initial_guess to set the initial guess for the nonlinear solver
                 *
                 *  Parameter \c coefficientName has default value set to "default", which is only used when setting
                 *  initial guess and allows us to call the function as
                 *  \code
                 *  setCoefficient ("initial_guess", function);
                 *  \endcode
                 *  without specifying the coefficient name (which would not be used anyway).
                 *  Note that when setting the initial guess, even though parameter \c coefficientValue is of type
                 *  \c boost::shared_ptr<const dolfin::GenericFunction>, the real type of the object pointed by 
                 *  \c coefficientValue can only be \c dolfin::Expression or \c dolfin::Function. This is because
                 *  \c setCoefficient will call the assignement operator of class \c dolfin::Function, which only accepts
                 *  the two types mentioned before as input arguments.
                 *  See \c AbstractDifferentialProblem documentation for more details on the function
                 */
                virtual void setCoefficient (const std::string& coefficientType, 
                                             const boost::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                             const std::string& coefficientName);

                //! Set coefficient [2]. Override of virtual function in \c AbstractDifferentialProblem.
                /*!
                 *  Possible values for \c coefficientType are:
                 *  \li residual_form to set the coefficient in the residual form
                 *  \li jacobian_form to set the coefficient in the jacobian form
                 *  
                 *  See \c AbstractDifferentialProblem documentation for more details on the function
                 */
                virtual void setCoefficient (const std::string& coefficientType,
                                             const boost::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                             const std::size_t& coefficientNumber);

                //! Set integration subdomains for the forms. Override of virtual function in \c AbstractDifferentialProblem
                /*! 
                 *  Possible values for \c coefficientType are:
                 *  \li residual_form to set the integration subdomain in the residual form
                 *  \li jacobian_form to set the integration subdomain in the jacobian form
                 *  
                 *  See \c AbstractDifferentialProblem documentation for more details on the function
                 */
                virtual void setIntegrationSubdomains (const std::string& formType,
                                                       boost::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                                       const control_problem::SubdomainType& subdomainType);
                
                /******************* METHODS *******************/

                //! Solve problem
                /*!
                 * This method solves the problem defined. It uses the private members' value to set the problem and then
                 * stores the solution in the private member \c solution_
                 */
                virtual void solve ();

                //! Solve problem specifying flag
                /*!
                 *  \param solverParameters object of type dolfin::parameters that contain the parameters to
                 *  be used for the non linear solver
                 */
                void solve (const dolfin::Parameters& solverParameters);

                //! Clone method. Overrides method in \c AbstractDifferentialProblem
                virtual control_problem::NonlinearDifferentialProblem <T_ResidualForm, T_JacobianForm>*
                    clone () const;
                
                // ---------------------------------------------------------------------------------------------//

            protected:
                //! The residual form
                T_ResidualForm residualForm_;

                //! The jacobian form
                T_JacobianForm jacobianForm_;
                
                // ---------------------------------------------------------------------------------------------//

            private:
        };



    // ==============================================================================================//
    // ==================================== IMPLEMENTATION ==========================================//
    // ==============================================================================================//


    /******************* CONSTRUCTORS *******************/

    template <class T_ResidualForm, class T_JacobianForm>
        NonlinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
        NonlinearDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
                                      const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
                                      const std::string& residualFormSolutionName,
                                      const std::string& jacobianFormSolutionName = "") : 
            AbstractDifferentialProblem (mesh, functionSpace),
            residualForm_ (*functionSpace),
            jacobianForm_ (*functionSpace, *functionSpace)
        { 
            dolfin::begin (dolfin::DBG, "Building NonlinearDifferentialProblem...");
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "nonlinear");
            parameters.add (dolfin::NonlinearVariationalSolver::default_parameters ());
            parameters.add ("residual_form_solution_name", residualFormSolutionName);
            if (jacobianFormSolutionName.empty ())
            {
                parameters.add ("jacobian_form_solution_name", residualFormSolutionName);
            }
            else
            {
                parameters.add ("jacobian_form_solution_name", jacobianFormSolutionName);
            }
            parameters.add ("solver_parameters_set_name", "nonlinear_variational_solver");
            
            dolfin::log (dolfin::DBG, "Setting initial guess...");
            boost::shared_ptr <dolfin::Function> tmpSolution (new dolfin::Function (solution_));
            if (jacobianFormSolutionName.empty ())
            {
                residualForm_.set_coefficient (residualFormSolutionName, tmpSolution);
                jacobianForm_.set_coefficient (residualFormSolutionName, tmpSolution);
            }
            else
            {
                residualForm_.set_coefficient (residualFormSolutionName, tmpSolution);
                jacobianForm_.set_coefficient (jacobianFormSolutionName, tmpSolution);
            }
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "NonlinearDifferentialProblem created");
        }



    template <class T_ResidualForm, class T_JacobianForm>
        NonlinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
        NonlinearDifferentialProblem (const dolfin::Mesh& mesh, 
                                      const dolfin::FunctionSpace& functionSpace,
                                      const std::string& residualFormSolutionName,
                                      const std::string& jacobianFormSolutionName = "") : 
            AbstractDifferentialProblem (mesh, functionSpace),
            residualForm_ (functionSpace),
            jacobianForm_ (functionSpace, functionSpace)
        { 
            dolfin::begin (dolfin::DBG, "Building NonlinearDifferentialProblem...");
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "nonlinear");
            parameters.add (dolfin::NonlinearVariationalSolver::default_parameters ());
            parameters.add ("residual_form_solution_name", residualFormSolutionName);
            if (jacobianFormSolutionName.empty ())
            {
                parameters.add ("jacobian_form_solution_name", residualFormSolutionName);
            }
            else
            {
                parameters.add ("jacobian_form_solution_name", jacobianFormSolutionName);
            }
            parameters.add ("solver_parameters_set_name", "nonlinear_variational_solver");
            
            dolfin::log (dolfin::DBG, "Setting initial guess...");
            boost::shared_ptr <dolfin::Function> tmpSolution (new dolfin::Function (solution_));
            if (jacobianFormSolutionName.empty ())
            {
                residualForm_.set_coefficient (residualFormSolutionName, tmpSolution);
                jacobianForm_.set_coefficient (residualFormSolutionName, tmpSolution);
            }
            else
            {
                residualForm_.set_coefficient (residualFormSolutionName, tmpSolution);
                jacobianForm_.set_coefficient (jacobianFormSolutionName, tmpSolution);
            }
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "NonlinearDifferentialProblem created");
        }
            

    
    template <class T_ResidualForm, class T_JacobianForm>
        NonlinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
        NonlinearDifferentialProblem (dolfin::Mesh&& mesh, 
                                      dolfin::FunctionSpace&& functionSpace,
                                      const std::string& residualFormSolutionName,
                                      const std::string& jacobianFormSolutionName = "") : 
            AbstractDifferentialProblem (mesh, functionSpace),
            residualForm_ (functionSpace),
            jacobianForm_ (functionSpace, functionSpace)
        { 
            dolfin::begin (dolfin::DBG, "Building NonlinearDifferentialProblem...");
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "nonlinear");
            parameters.add (dolfin::NonlinearVariationalSolver::default_parameters ());
            parameters.add ("residual_form_solution_name", residualFormSolutionName);
            if (jacobianFormSolutionName.empty ())
            {
                parameters.add ("jacobian_form_solution_name", residualFormSolutionName);
            }
            else
            {
                parameters.add ("jacobian_form_solution_name", jacobianFormSolutionName);
            }
            parameters.add ("solver_parameters_set_name", "nonlinear_variational_solver");
            
            dolfin::log (dolfin::DBG, "Setting initial guess...");
            boost::shared_ptr <dolfin::Function> tmpSolution (new dolfin::Function (solution_));
            if (jacobianFormSolutionName.empty ())
            {
                residualForm_.set_coefficient (residualFormSolutionName, tmpSolution);
                jacobianForm_.set_coefficient (residualFormSolutionName, tmpSolution);
            }
            else
            {
                residualForm_.set_coefficient (residualFormSolutionName, tmpSolution);
                jacobianForm_.set_coefficient (jacobianFormSolutionName, tmpSolution);
            }
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "NonlinearDifferentialProblem created");
        }
            


    template <class T_ResidualForm, class T_JacobianForm>
        NonlinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
        NonlinearDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
                                      const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
                                      const T_ResidualForm& residualForm,
                                      const T_JacobianForm& jacobianForm,
                                      const std::string& residualFormSolutionName,
                                      const std::string& jacobianFormSolutionName = "") : 
            AbstractDifferentialProblem (mesh, functionSpace),
            residualForm_ (residualForm),
            jacobianForm_ (jacobianForm)
        { 
            dolfin::begin (dolfin::DBG, "Building NonlinearDifferentialProblem...");
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "nonlinear");
            parameters.add (dolfin::NonlinearVariationalSolver::default_parameters ());
            parameters.add ("residual_form_solution_name", residualFormSolutionName);
            if (jacobianFormSolutionName.empty ())
            {
                parameters.add ("jacobian_form_solution_name", residualFormSolutionName);
            }
            else
            {
                parameters.add ("jacobian_form_solution_name", jacobianFormSolutionName);
            }
            parameters.add ("solver_parameters_set_name", "nonlinear_variational_solver");
            
            dolfin::log (dolfin::DBG, "Setting initial guess...");
            boost::shared_ptr <dolfin::Function> tmpSolution (new dolfin::Function (solution_));
            if (jacobianFormSolutionName.empty ())
            {
                residualForm_.set_coefficient (residualFormSolutionName, tmpSolution);
                jacobianForm_.set_coefficient (residualFormSolutionName, tmpSolution);
            }
            else
            {
                residualForm_.set_coefficient (residualFormSolutionName, tmpSolution);
                jacobianForm_.set_coefficient (jacobianFormSolutionName, tmpSolution);
            }
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "NonlinearDifferentialProblem created");
        }

    

    template <class T_ResidualForm, class T_JacobianForm>
        NonlinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
        NonlinearDifferentialProblem (const dolfin::Mesh& mesh, 
                                      const dolfin::FunctionSpace& functionSpace,
                                      const T_ResidualForm& residualForm,
                                      const T_JacobianForm& jacobianForm,
                                      const std::string& residualFormSolutionName,
                                      const std::string& jacobianFormSolutionName = "") : 
            AbstractDifferentialProblem (mesh, functionSpace),
            residualForm_ (residualForm),
            jacobianForm_ (jacobianForm)
        { 
            dolfin::begin (dolfin::DBG, "Building NonlinearDifferentialProblem...");
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "nonlinear");
            parameters.add (dolfin::NonlinearVariationalSolver::default_parameters ());
            parameters.add ("residual_form_solution_name", residualFormSolutionName);
            if (jacobianFormSolutionName.empty ())
            {
                parameters.add ("jacobian_form_solution_name", residualFormSolutionName);
            }
            else
            {
                parameters.add ("jacobian_form_solution_name", jacobianFormSolutionName);
            }
            parameters.add ("solver_parameters_set_name", "nonlinear_variational_solver");
            
            dolfin::log (dolfin::DBG, "Setting initial guess...");
            boost::shared_ptr <dolfin::Function> tmpSolution (new dolfin::Function (solution_));
            if (jacobianFormSolutionName.empty ())
            {
                residualForm_.set_coefficient (residualFormSolutionName, tmpSolution);
                jacobianForm_.set_coefficient (residualFormSolutionName, tmpSolution);
            }
            else
            {
                residualForm_.set_coefficient (residualFormSolutionName, tmpSolution);
                jacobianForm_.set_coefficient (jacobianFormSolutionName, tmpSolution);
            }
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "NonlinearDifferentialProblem created");
        }
            


    template <class T_ResidualForm, class T_JacobianForm>
        NonlinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
        NonlinearDifferentialProblem (dolfin::Mesh&& mesh, 
                                      dolfin::FunctionSpace&& functionSpace,
                                      T_ResidualForm&& residualForm,
                                      T_JacobianForm&& jacobianForm,
                                      const std::string& residualFormSolutionName,
                                      const std::string& jacobianFormSolutionName = "") : 
            AbstractDifferentialProblem (mesh, functionSpace),
            residualForm_ (residualForm),
            jacobianForm_ (jacobianForm)
        { 
            dolfin::begin (dolfin::DBG, "Building NonlinearDifferentialProblem...");
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "nonlinear");
            parameters.add (dolfin::NonlinearVariationalSolver::default_parameters ());
            parameters.add ("residual_form_solution_name", residualFormSolutionName);
            if (jacobianFormSolutionName.empty ())
            {
                parameters.add ("jacobian_form_solution_name", residualFormSolutionName);
            }
            else
            {
                parameters.add ("jacobian_form_solution_name", jacobianFormSolutionName);
            }
            parameters.add ("solver_parameters_set_name", "nonlinear_variational_solver");
            
            dolfin::log (dolfin::DBG, "Setting initial guess...");
            
            boost::shared_ptr <dolfin::Function> tmpSolution (new dolfin::Function (solution_));
            if (jacobianFormSolutionName.empty ())
            {
                residualForm_.set_coefficient (residualFormSolutionName, tmpSolution);
                jacobianForm_.set_coefficient (residualFormSolutionName, tmpSolution);
            }
            else
            {
                residualForm_.set_coefficient (residualFormSolutionName, tmpSolution);
                jacobianForm_.set_coefficient (jacobianFormSolutionName, tmpSolution);
            }
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "NonlinearDifferentialProblem created");
        }
    



    /******************* DESTRUCTOR *******************/

    // this is done for compatibility with gcc-4.6, which doesn't allow virtual members to be defualted in class body
    template <class T_ResidualForm, class T_JacobianForm>
        NonlinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
        ~NonlinearDifferentialProblem () = default;
    
    


    /******************* GETTERS *******************/

    template <class T_ResidualForm, class T_JacobianForm>
        const T_ResidualForm& NonlinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
        residualForm () const
        {
            return residualForm_;
        }



    template <class T_ResidualForm, class T_JacobianForm>
        const T_JacobianForm& NonlinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
        jacobianForm () const
        {
            return jacobianForm_;
        }

    

    /******************* SETTERS *******************/

    template <class T_ResidualForm, class T_JacobianForm>
        void NonlinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
        setCoefficient (const std::string& coefficientType, 
                        const boost::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                        const std::string& coefficientName = "default")
        {
            if (coefficientType == "residual_form")
            {
                dolfin::log (dolfin::DBG, "Setting residual form coefficient \"%s\"...", coefficientName.c_str ());
                residualForm_.set_coefficient (coefficientName, coefficientValue);
            }
            else if (coefficientType == "jacobian_form")
            {
                dolfin::log (dolfin::DBG, "Setting jacobian form coefficient \"%s\"...", coefficientName.c_str ());
                jacobianForm_.set_coefficient (coefficientName, coefficientValue);
            }
            else if (coefficientType == "initial_guess")
            {
                dolfin::log (dolfin::DBG, "Setting initial guess...");
                // check whether coefficientValue is a pointer to dolfin::Function or dolfin::Expression
                if (boost::dynamic_pointer_cast<const dolfin::Function> (coefficientValue) != nullptr)
                {
                    solution_ = *(boost::dynamic_pointer_cast<const dolfin::Function> (coefficientValue));
                }
                else if (boost::dynamic_pointer_cast<const dolfin::Expression> (coefficientValue) != nullptr)
                {
                    solution_ = *(boost::dynamic_pointer_cast<const dolfin::Expression> (coefficientValue));
                }
                else
                {
                    dolfin::warning ("Cannot set initial guess in nonlinear differential problem. Input argument is \
                                     neither a dolfin::Function nor a dolfin::Expression");
                }
            }
            else
            {
                dolfin::warning ("Cannot set coefficient in non linear differential problem. Coefficient type \"%s\" unknown",
                                 coefficientType.c_str ());
            }
        }



    template <class T_ResidualForm, class T_JacobianForm>
        void NonlinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
        setCoefficient (const std::string& coefficientType,
                        const boost::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                        const std::size_t& coefficientNumber)
        {
            if (coefficientType == "residual_form")
            {
                dolfin::log (dolfin::DBG, "Setting residual form coefficient number %d...", coefficientNumber);
                residualForm_.set_coefficient (coefficientNumber, coefficientValue);
            }
            else if (coefficientType == "jacobian_form")
            {
                dolfin::log (dolfin::DBG, "Setting jacobian form coefficient number %d...", coefficientNumber);
                jacobianForm_.set_coefficient (coefficientNumber, coefficientValue);
            }
            else
            {
                dolfin::warning ("Cannot set coefficient in non linear differential problem. Coefficient type \"%s\" unknown",
                                 coefficientType.c_str ());
            }
        }



    template <class T_ResidualForm, class T_JacobianForm>
        void NonlinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
        setIntegrationSubdomains (const std::string& formType,
                                  boost::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                  const control_problem::SubdomainType& subdomainType)
        {
            if (formType == "residual_form")
            {
                if (subdomainType == control_problem::SubdomainType::INTERNAL_CELLS)
                {
                    dolfin::log (dolfin::DBG, "Setting residual form integration subdomain on INTERNAL_CELLS...");
                    residualForm_.set_cell_domains (meshFunction);
                }
                else if (subdomainType == control_problem::SubdomainType::INTERNAL_FACETS)
                {
                    dolfin::log (dolfin::DBG, "Setting residual form integration subdomain on INTERNAL_FACETS...");
                    residualForm_.set_interior_facet_domains (meshFunction);
                }
                else if (subdomainType == control_problem::SubdomainType::BOUNDARY_FACETS)
                {
                    dolfin::log (dolfin::DBG, "Setting residual form integration subdomain on BOUNDARY_FACETS...");
                    residualForm_.set_exterior_facet_domains (meshFunction);
                }
                else
                {
                    dolfin::warning ("unknown subdomain type requested while trying to apply mesh function to residual form"); 
                }
            }
            else if (formType == "jacobian_form")
            {
                if (subdomainType == control_problem::SubdomainType::INTERNAL_CELLS)
                {
                    dolfin::log (dolfin::DBG, "Setting jacobian form integration subdomain on INTERNAL_CELLS...");
                    jacobianForm_.set_cell_domains (meshFunction);
                }
                else if (subdomainType == control_problem::SubdomainType::INTERNAL_FACETS)
                {
                    dolfin::log (dolfin::DBG, "Setting jacobian form integration subdomain on INTERNAL_FACETS...");
                    jacobianForm_.set_interior_facet_domains (meshFunction);
                }
                else if (subdomainType == control_problem::SubdomainType::BOUNDARY_FACETS)
                {
                    dolfin::log (dolfin::DBG, "Setting jacobian form integration subdomain on BOUNDARY_FACETS...");
                    jacobianForm_.set_exterior_facet_domains (meshFunction);
                }
                else
                {
                    dolfin::warning ("unknown subdomain type requested while trying to apply mesh function to jacobian form"); 
                }
            }
            else
            {
                dolfin::warning ("Cannot set integration subdomain in linear differential problem. Form type \"%s\" unknown",
                                 formType.c_str ());
            }

        }



    template <class T_ResidualForm, class T_JacobianForm>
        void NonlinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
        solve () 
        {
            dolfin::begin (dolfin::DBG, "Solving problem...");
            
            dolfin::log (dolfin::DBG, "Creating temporary vectors of dolfin::DirichletBC pointers...");
            // create vector of POINTERS to DirichletBC. This is needed to call the function dolfin::solve
            std::vector<const dolfin::DirichletBC*> tmpDirichletBCs (dirichletBCs_.size (), nullptr);
            for (std::size_t i = 0; i < dirichletBCs_.size (); ++i)
            {
                tmpDirichletBCs[i] = &dirichletBCs_[i];
            }
            // note that when tmpDirichletBCs gets destroyd (upon exit of the function) the vector distructor will
            // call DirichletBC* destructor, which means that the pointers will be destroyed but the object they point
            // to will not
            
            // now solve non linear problem and store solution in solution_.
            // We need to check whether tmpDirichletBCs is empty, and call the appropriate function
            // Parameters passed to the solver are those stored in private member "parameters" in the set identified
            // by parameter "solver_parameters_set_name"
            std::string solverParametersSetName = parameters ["solver_parameters_set_name"];
            dolfin::log (dolfin::DBG, "Solver parameters set name is: %s", solverParametersSetName.c_str ());
            if (tmpDirichletBCs.size () != 0)
            {
                if (dolfin::get_log_level () > dolfin::DBG)
                {
                    dolfin::end ();
                }
                dolfin::solve (residualForm_ == 0, solution_, tmpDirichletBCs, jacobianForm_, 
                               parameters (solverParametersSetName));
            }
            else
            {
                if (dolfin::get_log_level () > dolfin::DBG)
                {
                    dolfin::end ();
                }
                dolfin::solve (residualForm_ == 0, solution_, jacobianForm_, 
                               parameters (solverParametersSetName));
            }
            
            if (dolfin::get_log_level () <= dolfin::DBG)
            {
                dolfin::end ();
            }
        }



    template <class T_ResidualForm, class T_JacobianForm>
        void NonlinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
        solve (const dolfin::Parameters& solverParameters)
        {
            dolfin::begin (dolfin::DBG, "Solving problem...");
            
            dolfin::log (dolfin::DBG, "Creating temporary vectors of dolfin::DirichletBC pointers...");
            // create vector of POINTERS to DirichletBC. This is needed to call the function dolfin::solve
            std::vector<const dolfin::DirichletBC*> tmpDirichletBCs (dirichletBCs_.size (), nullptr);
            for (std::size_t i = 0; i < dirichletBCs_.size (); ++i)
            {
                tmpDirichletBCs[i] = &dirichletBCs_[i];
            }
            // note that when tmpDirichletBCs gets destroyd (upon exit of the function) the vector distructor will
            // call DirichletBC* destructor, which means that the pointers will be destroyed but the object they point
            // to will not
            
            // now solve non linear problem and store solution in solution_.
            // We need to check whether tmpDirichletBCs is empty, and call the appropriate function
            // Parameters passed to the solver are those stored in the input variable
            if (tmpDirichletBCs.size () != 0)
            {
                if (dolfin::get_log_level () > dolfin::DBG)
                {
                    dolfin::end ();
                }
                dolfin::solve (residualForm_ == 0, solution_, tmpDirichletBCs, jacobianForm_, solverParameters);
            }
            else
            {
                if (dolfin::get_log_level () > dolfin::DBG)
                {
                    dolfin::end ();
                }
                dolfin::solve (residualForm_ == 0, solution_, jacobianForm_, solverParameters);
            }
            
            if (dolfin::get_log_level () <= dolfin::DBG)
            {
                dolfin::end ();
            }
        }
    
    
    
    template <class T_ResidualForm, class T_JacobianForm>
        control_problem::NonlinearDifferentialProblem <T_ResidualForm, T_JacobianForm>*
        NonlinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
        clone () const
        {
            dolfin::begin (dolfin::DBG, "Cloning object...");
            
            if (dolfin::get_log_level () > dolfin::DBG)
            {
                dolfin::end ();
            }
            
            dolfin::log (dolfin::DBG, "Creating new object of type NonlinearDifferentialProblem...");
            // create new object
            control_problem::NonlinearDifferentialProblem <T_ResidualForm, T_JacobianForm>* clonedProblem 
                (new control_problem::NonlinearDifferentialProblem <T_ResidualForm, T_JacobianForm> 
                        ( this->mesh_,
                          this->functionSpace_,
                          this->residualForm_, 
                          this->jacobianForm_,
                         (this->parameters) ["residual_form_solution_name"],
                         (this->parameters) ["jacobian_form_solution_name"]
                        )
                );
            
            //copy dirichlet boundary conditions
            dolfin::log (dolfin::DBG, "Copying Dirichlet boundary conditions...");
            for (auto i : this->dirichletBCs_)
            {
                clonedProblem->addDirichletBC (i);
            }
            
            // clear parameters set of newly created object so that it can be populated by the parameters of the object
            // being created
            dolfin::log (dolfin::DBG, "Copying parameters to new object...");
            clonedProblem->parameters.clear ();
            clonedProblem->parameters = this->parameters;
            
            // copy solution
            dolfin::log (dolfin::DBG, "Copying solution...");
            clonedProblem->solution_ = this->solution_;
            
            if (dolfin::get_log_level () <= dolfin::DBG)
            {
                dolfin::end ();
            }
            
            return clonedProblem;
        }
}
#endif

