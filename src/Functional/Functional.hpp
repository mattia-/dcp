#ifndef SRC_FUNCTIONAL_FUNCTIONAL_HPP_INCLUDE_GUARD
#define SRC_FUNCTIONAL_FUNCTIONAL_HPP_INCLUDE_GUARD

#include <dolfin.h>
#include <string>
#include <memory>
#include <Utils/SubdomainType.hpp>

namespace controlproblem
{
    /*! \class Functional Functional.hpp
     *  \brief Class for functionals.
     *
     *  This class represents a generic functional and contains both its representation and its gradient.
     *  The class is template-ized over the type of the functional, which must be derived from \c dolfin::Form (and will
     *  typically defined in a ufl file with the keyword \c forms).
     *  
     *  Template arguments are:
     *  \arg T_FunctionalForm_ the functional form type
     */

    template <class T_FunctionalForm_>
        class Functional
        {
            // ---------------------------------------------------------------------------------------------//  

            public:
                /******************* TYPEDEFS *******************/
                typedef T_FunctionalForm_ T_FunctionalForm;
                
                
                /******************* CONSTRUCTORS *******************/
                //! Default constructor is deleted. The class is not default constructable.
                Functional () = delete;

                //!  Constructor with shared pointers [1]
                /*!
                 *  \param mesh the mesh over which the functional is defined as a \c const \c std::shared_ptr to 
                 *  \c dolfin::Mesh
                 *  The stored mesh's ownership will be shared between the object and the input argument.
                 *  The functional form  will be created too, calling the constructor which takes the mesh
                 *  as input.
                 */
                Functional (const std::shared_ptr<dolfin::Mesh> mesh);
                

                //! Constructor with references [1]
                /*!
                 *  \param mesh the mesh as a \c const \c dolfin::Mesh&
                 *  The stored mesh's ownership will be unique to the object, since the protected member is 
                 *  initialized using the \c new operator and mesh's copy constructor.
                 *  The functional form will be created too, calling the constructor which takes the mesh
                 *  as input.
                 */
                Functional (const dolfin::Mesh& mesh);

                //! Constructor with rvalue references [1]
                /*!
                 *  \param mesh the mesh as a \c dolfin::Mesh&&
                 *  The stored mesh's ownership will be unique to the object, since the protected member is 
                 *  initialized using the \c new operator and mesh's move constructor
                 *  The functional form will be created too, calling the constructor which takes the mesh
                 *  as input.
                 */
                Functional (dolfin::Mesh&& mesh);

                
                //!  Constructor with shared pointers [2]
                /*!
                 *  \param mesh the mesh as a \c const \c std::shared_ptr to \c dolfin::Mesh
                 *  The stored mesh's ownership will be shared between the object and the input argument.
                 *  The functional form will be created too, calling the constructor which takes the mesh
                 *  as input.
                 */
                Functional (const std::shared_ptr<dolfin::Mesh> mesh);

                //! Constructor with references [2]
                /*!
                 *  \param mesh the mesh as a \c const \c dolfin::Mesh&
                 *  The stored mesh's ownership will be unique to the object, since the protected member is 
                 *  initialized using the \c new operator and mesh's copy constructor
                 *  The functional form will be created too, calling the constructor which takes the mesh
                 *  as input.
                 */
                Functional (const dolfin::Mesh& mesh);

                //! Constructor with rvalue references [3]
                /*!
                 *  \param mesh the mesh as a \c dolfin::Mesh&&
                 *  The stored mesh's ownership will be unique to the object, since the protected member is 
                 *  initialized using the \c new operator and mesh's move constructor
                 *  The functional form will be created too, calling the constructor which takes the mesh
                 *  as input.
                 */
                Functional (dolfin::Mesh&& mesh);

                /******************* DESTRUCTOR *******************/
                
                //! Destructor
                /*! 
                 *  Default destructor, since members of the class are trivially 
                 *  destructible.
                 */
                 ~Functional () = default;

                
                /******************* GETTERS *******************/
                //! Get const reference to the functional form
                /*! 
                 *  \return a const reference to the functional form
                 */
                const T_FunctionalForm& functionalForm () const;


                /******************* SETTERS *******************/
                
                //! Set coefficient [1]
                /*!
                 *  Possible values for \c coefficientType are:
                 *  \li bilinear_form to set the coefficient in the bilinear form
                 *  \li linear_form to set the coefficient in the linear form
                 *  
                 *  See \c AbstractDifferentialProblem documentation for more details on the function
                 */
                 void setCoefficient (const std::string& coefficientType, 
                                      const boost::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                      const std::string& coefficientName);

                //! Set coefficient [2]. Override of  function in \c AbstractDifferentialProblem.
                /*!
                 *  Possible values for \c coefficientType are:
                 *  \li bilinear_form to set the coefficient in the bilinear form
                 *  \li linear_form to set the coefficient in the linear form
                 *  
                 *  See \c AbstractDifferentialProblem documentation for more details on the function
                 */
                 void setCoefficient (const std::string& coefficientType,
                                      const boost::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                      const std::size_t& coefficientNumber);

                 //! Set integration subdomains for the forms. Override of  function in \c AbstractDifferentialProblem
                 /*! 
                  *  Possible values for \c coefficientType are:
                  *  \li bilinear_form to set the integration subdomain in the bilinear form
                  *  \li linear_form to set the integration subdomain in the linear form
                  *  
                  *  See \c AbstractDifferentialProblem documentation for more details on the function
                  */
                 void setIntegrationSubdomains (const std::string& formType,
                                                boost::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                                const controlproblem::SubdomainType& subdomainType);

                 //! Add Dirichlet boundary condition to the problem [1]. Overrides method in \c AbstractDifferentialProblem
                 /*!
                  *  This method adds to the base class method the setting of parameter \c system_is_assembled to false.
                  *  \param dirichletCondition a const reference to the dirichlet boundary condition to be added to the problem
                  */
                 void addDirichletBC (const dolfin::DirichletBC& dirichletCondition);

                 //! Add Dirichlet boundary condition to the problem [2]. Overrides method in \c AbstractDifferentialProblem
                 /*!
                  *  This method adds to the base class method the setting of parameter \c system_is_assembled to false.
                  *  \param dirichletCondition a rvalue reference to the dirichlet boundary condition to be added to the problem
                  */
                 void addDirichletBC (dolfin::DirichletBC&& dirichletCondition);

                 //! Remove Dirichlet boundary condition with given position. Overrides method in \c AbstractDifferentialProblem
                 /*!
                  *  This method adds to the base class method the setting of parameter \c system_is_assembled to false.
                  *  \param i the position in the vector of the boundary condition to be removed.
                  *            If i is greater than the size of the vector, nothing is removed.
                  */
                 void removeDirichletBC (const std::vector<dolfin::DirichletBC>::iterator& i);

                 //! Method to update class members. It checks for differences between desired and current solver parameters
                 //! and creates a new solver, setting also the proper parameters
                 void update ();


                 /******************* METHODS *******************/

                 //! Solve problem
                 /*!
                  *  This method solves the problem defined. It uses the private members' value to set the problem and then
                  *  stores the solution in the private member \c solution_. Note that it checks the protected member
                  *  \c parameters to decided whether problem's matrix and vector should be reassembled and if
                  *  the values of the parameters "desired_solver_type", "desired_solver_method" and 
                  *  "desired_solver_preconditioner" match the values of "current_solver_type", "current_solver_method"
                  *  and "current_solver_preconditioner". If they differ, it calls \c createSolver()
                  */
                 void solve ();

                 //! Solve problem specifying flag
                 /*!
                  *  \param mustReassemble true if the system operators (matrix and right hand side vector)
                  *         should be reassembled. It is false by default.
                  */
                void solve (const bool& mustReassemble);

                //! Clone method. Overrides method in \c AbstractDifferentialProblem
                 controlproblem::Functional<T_FunctionalForm>*
                    clone () const;

                // ---------------------------------------------------------------------------------------------//

            protected:
                //! Creates a linear solver object of type passed as input. 
                /*!
                 *  The solver will be created using the parameters "desired_solver_type", "desired_solver_method"
                 *  and "desired_solver_preconditioner" set in the protected member \c parameters. It will
                 *  also set the parameters "current_solver_type", "current_solver_method" and 
                 *  "current_solver_preconditioner" in the same set of parameters substituting the current values with 
                 *  the values being used to create the solver
                 */
                std::unique_ptr<dolfin::GenericLinearSolver> createSolver ();
                
                //! The bilinear form
                T_BilinearForm bilinearForm_;

                //! The linear form
                T_LinearForm linearForm_;

                //! The solver
                /*! 
                 *  We use a pointer so that polymorphism can be applied.
                 */
                std::unique_ptr<dolfin::GenericLinearSolver> solver_;
                
                //! Matrix to hold the problem's discrete operator
                /*!
                 *  We use a boost::shared_ptr for compatibility with FEniCS, version 1.3.0.
                 *  Note that there is no write access to this variable from outside the class,
                 *  so it is guaranteed to be a unique shared_ptr
                 */
                boost::shared_ptr<dolfin::Matrix> problemMatrix_;

                //! Vector to store the right hand side of the discrete problem
                dolfin::Vector rhsVector_;

                // ---------------------------------------------------------------------------------------------//

            private:
        };



    // ============================================================================================== //
    // ==================================== IMPLEMENTATION ========================================== //
    // ============================================================================================== //


    /******************* CONSTRUCTORS *******************/

    template <class T_FunctionalForm>
        Functional<T_FunctionalForm>::
        Functional (const std::shared_ptr<dolfin::Mesh> mesh, 
                                   const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
                                   const std::string& solverType = "lu_solver",
                                   const std::string& solverMethod = "default",
                                   const std::string& solverPreconditioner = "default") :
            AbstractDifferentialProblem (mesh, functionSpace),
            bilinearForm_ (*functionSpace, *functionSpace),
            linearForm_ (*functionSpace),
            solver_ (nullptr),
            problemMatrix_ (new dolfin::Matrix),
            rhsVector_ ()
        { 
            dolfin::begin (dolfin::DBG, "Building Functional...");
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "linear");
            parameters.add ("current_solver_type", solverType);
            parameters.add ("current_solver_method", solverMethod);
            parameters.add ("current_solver_preconditioner", solverPreconditioner);
            parameters.add ("desired_solver_type", solverType);
            parameters.add ("desired_solver_method", solverMethod);
            parameters.add ("desired_solver_preconditioner", solverPreconditioner);
            parameters.add ("system_is_assembled", false);
            parameters.add ("force_reassemble_system", false);  
            
            dolfin::begin (dolfin::DBG, "Creating solver...");
            solver_ = createSolver ();
            dolfin::end ();
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "Functional created");
        }



    template <class T_FunctionalForm>
        Functional<T_FunctionalForm>::
        Functional (const dolfin::Mesh& mesh, 
                                   const dolfin::FunctionSpace& functionSpace,
                                   const std::string& solverType = "lu_solver",
                                   const std::string& solverMethod = "default",
                                   const std::string& solverPreconditioner = "default") :
            AbstractDifferentialProblem (mesh, functionSpace),
            bilinearForm_ (functionSpace, functionSpace),
            linearForm_ (functionSpace),
            solver_ (nullptr),
            problemMatrix_ (new dolfin::Matrix),
            rhsVector_ ()
        { 
            dolfin::begin (dolfin::DBG, "Building Functional...");
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "linear");
            parameters.add ("current_solver_type", solverType);
            parameters.add ("current_solver_method", solverMethod);
            parameters.add ("current_solver_preconditioner", solverPreconditioner);
            parameters.add ("desired_solver_type", solverType);
            parameters.add ("desired_solver_method", solverMethod);
            parameters.add ("desired_solver_preconditioner", solverPreconditioner);
            parameters.add ("system_is_assembled", false);
            parameters.add ("force_reassemble_system", false);  
            
            dolfin::begin (dolfin::DBG, "Creating solver...");
            solver_ = createSolver ();
            dolfin::end ();
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "Functional created");
        }



    template <class T_FunctionalForm>
        Functional<T_FunctionalForm>::
        Functional (dolfin::Mesh&& mesh, 
                                   dolfin::FunctionSpace&& functionSpace,
                                   const std::string& solverType = "lu_solver",
                                   const std::string& solverMethod = "default",
                                   const std::string& solverPreconditioner = "default") :
            AbstractDifferentialProblem (mesh, functionSpace),
            bilinearForm_ (functionSpace, functionSpace),
            linearForm_ (functionSpace),
            solver_ (nullptr),
            problemMatrix_ (new dolfin::Matrix),
            rhsVector_ ()
        { 
            dolfin::begin (dolfin::DBG, "Building Functional...");
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "linear");
            parameters.add ("current_solver_type", solverType);
            parameters.add ("current_solver_method", solverMethod);
            parameters.add ("current_solver_preconditioner", solverPreconditioner);
            parameters.add ("desired_solver_type", solverType);
            parameters.add ("desired_solver_method", solverMethod);
            parameters.add ("desired_solver_preconditioner", solverPreconditioner);
            parameters.add ("system_is_assembled", false);
            parameters.add ("force_reassemble_system", false);
        
            dolfin::begin (dolfin::DBG, "Creating solver...");
            solver_ = createSolver ();
            dolfin::end ();
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "Functional created");
        }



    template <class T_FunctionalForm>
        Functional<T_FunctionalForm>::
        Functional (const std::shared_ptr<dolfin::Mesh> mesh, 
                                   const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
                                   const T_BilinearForm& bilinearForm,
                                   const T_LinearForm& linearForm,
                                   const std::string& solverType = "lu_solver",
                                   const std::string& solverMethod = "default",
                                   const std::string& solverPreconditioner = "default") :
            AbstractDifferentialProblem (mesh, functionSpace),
            bilinearForm_ (bilinearForm),
            linearForm_ (linearForm),
            solver_ (nullptr),
            problemMatrix_ (new dolfin::Matrix),
            rhsVector_ ()
        { 
            dolfin::begin (dolfin::DBG, "Building Functional...");
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "linear");
            parameters.add ("current_solver_type", solverType);
            parameters.add ("current_solver_method", solverMethod);
            parameters.add ("current_solver_preconditioner", solverPreconditioner);
            parameters.add ("desired_solver_type", solverType);
            parameters.add ("desired_solver_method", solverMethod);
            parameters.add ("desired_solver_preconditioner", solverPreconditioner);
            parameters.add ("system_is_assembled", false);
            parameters.add ("force_reassemble_system", false);
            
            dolfin::begin (dolfin::DBG, "Creating solver...");
            solver_ = createSolver ();
            dolfin::end ();
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "Functional created");
        }



    template <class T_FunctionalForm>
        Functional<T_FunctionalForm>::
        Functional (const dolfin::Mesh& mesh, 
                                   const dolfin::FunctionSpace& functionSpace,
                                   const T_BilinearForm& bilinearForm,
                                   const T_LinearForm& linearForm,
                                   const std::string& solverType = "lu_solver",
                                   const std::string& solverMethod = "default",
                                   const std::string& solverPreconditioner = "default") :
            AbstractDifferentialProblem (mesh, functionSpace),
            bilinearForm_ (bilinearForm),
            linearForm_ (linearForm),
            solver_ (nullptr),
            problemMatrix_ (new dolfin::Matrix),
            rhsVector_ ()
        { 
            dolfin::begin (dolfin::DBG, "Building Functional...");
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "linear");
            parameters.add ("current_solver_type", solverType);
            parameters.add ("current_solver_method", solverMethod);
            parameters.add ("current_solver_preconditioner", solverPreconditioner);
            parameters.add ("desired_solver_type", solverType);
            parameters.add ("desired_solver_method", solverMethod);
            parameters.add ("desired_solver_preconditioner", solverPreconditioner);
            parameters.add ("system_is_assembled", false);
            parameters.add ("force_reassemble_system", false);  
            
            dolfin::begin (dolfin::DBG, "Creating solver...");
            solver_ = createSolver ();
            dolfin::end ();
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "Functional created");
        }



    template <class T_FunctionalForm>
        Functional<T_FunctionalForm>::
        Functional (dolfin::Mesh&& mesh, 
                                   dolfin::FunctionSpace&& functionSpace,
                                   T_BilinearForm&& bilinearForm,
                                   T_LinearForm&& linearForm,
                                   const std::string& solverType = "lu_solver",
                                   const std::string& solverMethod = "default",
                                   const std::string& solverPreconditioner = "default") :
            AbstractDifferentialProblem (mesh, functionSpace),
            bilinearForm_ (std::move (bilinearForm)),
            linearForm_ (std::move (linearForm)),
            solver_ (nullptr),
            problemMatrix_ (new dolfin::Matrix),
            rhsVector_ ()
        {
            dolfin::begin (dolfin::DBG, "Building Functional...");
            
            dolfin::log (dolfin::DBG, "Setting up parameters...");
            parameters.add ("problem_type", "linear");
            parameters.add ("current_solver_type", solverType);
            parameters.add ("current_solver_method", solverMethod);
            parameters.add ("current_solver_preconditioner", solverPreconditioner);
            parameters.add ("desired_solver_type", solverType);
            parameters.add ("desired_solver_method", solverMethod);
            parameters.add ("desired_solver_preconditioner", solverPreconditioner);
            parameters.add ("system_is_assembled", false);
            parameters.add ("force_reassemble_system", false);  
            
            dolfin::begin (dolfin::DBG, "Creating solver...");
            solver_ = createSolver ();
            dolfin::end ();
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "Functional created");
        }

    
    /******************* DESTRUCTOR *******************/

    // this is done for compatibility with gcc-4.6, which doesn't allow  members to be defualted in class body
    template <class T_FunctionalForm>
        Functional<T_FunctionalForm>::
        ~Functional () = default;
    
    


    /******************* GETTERS *******************/

    template <class T_FunctionalForm>
        const T_BilinearForm& Functional<T_FunctionalForm>::
        bilinearForm () const
        {
            return bilinearForm_;
        }



    template <class T_FunctionalForm>
        const T_LinearForm& Functional<T_FunctionalForm>::
        linearForm () const
        {
            return linearForm_;
        }



    template <class T_FunctionalForm>
        const dolfin::Matrix& Functional<T_FunctionalForm>::
        linearOperator () const
        {
            return *problemMatrix_;
        }



    template <class T_FunctionalForm>
        const dolfin::Vector& Functional<T_FunctionalForm>::
        rhs () const
        {
            return rhsVector_;
        }



    /******************* SETTERS *******************/

    template <class T_FunctionalForm>
        void Functional<T_FunctionalForm>::
        setCoefficient (const std::string& coefficientType, 
                        const boost::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                        const std::string& coefficientName)
        {
            if (coefficientType == "bilinear_form")
            {
                dolfin::log (dolfin::DBG, "Setting bilinear form coefficient \"%s\"...", coefficientName.c_str ());
                bilinearForm_.set_coefficient (coefficientName, coefficientValue);
                parameters ["system_is_assembled"] = false;
            }
            else if (coefficientType == "linear_form")
            {
                dolfin::log (dolfin::DBG, "Setting linear form coefficient \"%s\"...", coefficientName.c_str ());
                linearForm_.set_coefficient (coefficientName, coefficientValue);
                parameters ["system_is_assembled"] = false;
            }
            else
            {
                dolfin::warning ("Cannot set coefficient in linear differential problem. Form type \"%s\" unknown",
                                 coefficientType.c_str ());
            }

        }

    

    template <class T_FunctionalForm>
        void Functional<T_FunctionalForm>::
        setCoefficient (const std::string& coefficientType,
                        const boost::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                        const std::size_t& coefficientNumber)
        {
            if (coefficientType == "bilinear_form")
            {
                dolfin::log (dolfin::DBG, "Setting bilinear form coefficient number %d...", coefficientNumber);
                bilinearForm_.set_coefficient (coefficientNumber, coefficientValue);
                parameters ["system_is_assembled"] = false;
            }
            else if (coefficientType == "linear_form")
            {
                dolfin::log (dolfin::DBG, "Setting linear form coefficient number %d...", coefficientNumber);
                linearForm_.set_coefficient (coefficientNumber, coefficientValue);
                parameters ["system_is_assembled"] = false;
            }
            else
            {
                dolfin::warning ("Cannot set coefficient in linear differential problem. Form type \"%s\" unknown",
                                 coefficientType.c_str ());
            }
        }



    template <class T_FunctionalForm>
        void Functional<T_FunctionalForm>::
        setIntegrationSubdomains (const std::string& formType,
                                  boost::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                  const controlproblem::SubdomainType& subdomainType)
        {
            if (formType == "bilinear_form")
            {
                if (subdomainType == controlproblem::SubdomainType::INTERNAL_CELLS)
                {
                    dolfin::log (dolfin::DBG, "Setting bilinear form integration subdomain on INTERNAL_CELLS...");
                    bilinearForm_.set_cell_domains (meshFunction);
                    parameters ["system_is_assembled"] = false;
                }
                else if (subdomainType == controlproblem::SubdomainType::INTERNAL_FACETS)
                {
                    dolfin::log (dolfin::DBG, "Setting bilinear form integration subdomain on INTERNAL_FACETS...");
                    bilinearForm_.set_interior_facet_domains (meshFunction);
                    parameters ["system_is_assembled"] = false;
                }
                else if (subdomainType == controlproblem::SubdomainType::BOUNDARY_FACETS)
                {
                    dolfin::log (dolfin::DBG, "Setting bilinear form integration subdomain on BOUNDARY_FACETS...");
                    bilinearForm_.set_exterior_facet_domains (meshFunction);
                    parameters ["system_is_assembled"] = false;
                }
                else
                {
                    dolfin::warning ("unknown subdomain type requested while trying to apply mesh function to bilinear form"); 
                }
            }
            else if (formType == "linear_form")
            {
                if (subdomainType == controlproblem::SubdomainType::INTERNAL_CELLS)
                {
                    dolfin::log (dolfin::DBG, "Setting linear form integration subdomain on INTERNAL_CELLS...");
                    linearForm_.set_cell_domains (meshFunction);
                    parameters ["system_is_assembled"] = false;
                }
                else if (subdomainType == controlproblem::SubdomainType::INTERNAL_FACETS)
                {
                    dolfin::log (dolfin::DBG, "Setting linear form integration subdomain on INTERNAL_FACETS...");
                    linearForm_.set_interior_facet_domains (meshFunction);
                    parameters ["system_is_assembled"] = false;
                }
                else if (subdomainType == controlproblem::SubdomainType::BOUNDARY_FACETS)
                {
                    dolfin::log (dolfin::DBG, "Setting linear form integration subdomain on BOUNDARY_FACETS...");
                    linearForm_.set_exterior_facet_domains (meshFunction);
                }
                else
                {
                    dolfin::warning ("unknown subdomain type requested while trying to apply mesh function to linear form"); 
                }
            }
            else
            {
                dolfin::warning ("Cannot set integration subdomain in linear differential problem. Form type \"%s\" unknown",
                                 formType.c_str ());
            }

        }



    template <class T_FunctionalForm>
        void Functional<T_FunctionalForm>::
        addDirichletBC (const dolfin::DirichletBC& dirichletCondition)
        {
            dirichletBCs_.emplace_back (dirichletCondition);
            parameters ["system_is_assembled"] = false;
        }


    template <class T_FunctionalForm>
        void Functional<T_FunctionalForm>::
        addDirichletBC (dolfin::DirichletBC&& dirichletCondition)
        {
            dirichletBCs_.emplace_back (dirichletCondition);
            parameters ["system_is_assembled"] = false;
        }



    template <class T_FunctionalForm>
        void Functional<T_FunctionalForm>::
        removeDirichletBC (const std::vector<dolfin::DirichletBC>::iterator& i)
        {
            dirichletBCs_.erase (i);
            parameters ["system_is_assembled"] = false;
        }



    template <class T_FunctionalForm>
        void Functional<T_FunctionalForm>::
        update ()
        {
            // define auxiliary string variables
            std::string desiredSolverType = parameters ["desired_solver_type"];
            std::string desiredSolverMethod = parameters ["desired_solver_method"];
            std::string desiredSolverPreconditioner = parameters ["desired_solver_preconditioner"];
            
            std::string currentSolverType = parameters ["current_solver_type"];
            std::string currentSolverMethod = parameters ["current_solver_method"];
            std::string currentSolverPreconditioner = parameters ["current_solver_preconditioner"];
            
            bool needSolverUpdate = desiredSolverType != currentSolverType 
                                    || desiredSolverMethod != currentSolverMethod 
                                    || desiredSolverPreconditioner != currentSolverPreconditioner;
            
            if (needSolverUpdate)
            {
                dolfin::begin (dolfin::DBG, "Updating solver...");
                solver_ = createSolver ();
                dolfin::end ();
                
                parameters ["system_is_assembled"] = false;
            }
        }



    template <class T_FunctionalForm>
        void Functional<T_FunctionalForm>::
        solve () 
        {
            update ();
            
            // define auxiliary string variables
            bool systemIsAssembled = parameters ["system_is_assembled"];
            bool forceReassembleSystem = parameters ["force_reassemble_system"];
            bool needReassemble = !systemIsAssembled || forceReassembleSystem;
            
            if (needReassemble)
            {
                dolfin::begin (dolfin::DBG, "Assembling system...");
                
                dolfin::assemble (*problemMatrix_, bilinearForm_);
                dolfin::assemble (rhsVector_, linearForm_);
                for (auto i : dirichletBCs_)
                {
                    i.apply (*problemMatrix_, rhsVector_);
                }
                
                solver_ -> set_operator (problemMatrix_);
                
                parameters ["system_is_assembled"] = true;
                
                dolfin::end ();
            }
            
            dolfin::begin (dolfin::DBG, "Solving problem...");
            if (dolfin::get_log_level () > dolfin::DBG)
            {
                dolfin::end ();
            }
            
            solver_ -> solve (*solution_.vector (), rhsVector_);
            
            if (dolfin::get_log_level () <= dolfin::DBG)
            {
                dolfin::end ();
            }
        }



    template <class T_FunctionalForm>
        void Functional<T_FunctionalForm>::
        solve (const bool& mustReassemble)
        {
            // save current value of parameter force_reassemble_system
            bool forceReassembleSystemBackup = parameters ["force_reassemble_system"];
            
            // overwrite parameter force_reassemble_system with input value and solve problem
            parameters ["force_reassemble_system"] = mustReassemble;
            solve ();
            
            // restore value of parameter force_reassemble_system
            parameters ["force_reassemble_system"] = forceReassembleSystemBackup;
        }
    
    
    
    template <class T_FunctionalForm>
        controlproblem::Functional<T_FunctionalForm>*
        Functional<T_FunctionalForm>::
        clone () const
        {
            dolfin::begin (dolfin::DBG, "Cloning object...");
            
            if (dolfin::get_log_level () > dolfin::DBG)
            {
                dolfin::end ();
            }
            
            dolfin::log (dolfin::DBG, "Creating new object of type Functional...");
            
            // create new object
            controlproblem::Functional <T_FunctionalForm>*
                clonedProblem 
                (new controlproblem::Functional <T_FunctionalForm> 
                        (this->mesh_,
                         this->functionSpace_,
                         this->bilinearForm_, 
                         this->linearForm_
                        )
                );
            
            //copy dirichlet boundary conditions
            dolfin::log (dolfin::DBG, "Copying Dirichlet boundary conditions...");
            for (auto i : this->dirichletBCs_)
            {
                clonedProblem->addDirichletBC (i);
            }
            
            // clear parameters set of newly created object so that it can be populated by the parameters of the object
            // being created. Set "system_is_assembled" to false, though, because in order for it to work the newly
            // created problem will have to reassemble matrix and rhs vector
            dolfin::log (dolfin::DBG, "Copying parameters to new object...");
            clonedProblem->parameters.clear ();
            clonedProblem->parameters = this->parameters;
            clonedProblem->update ();
            clonedProblem->parameters ["system_is_assembled"] = false;
            
            // copy solution
            dolfin::log (dolfin::DBG, "Copying solution...");
            clonedProblem->solution_ = this->solution_;
            
            if (dolfin::get_log_level () <= dolfin::DBG)
            {
                dolfin::end ();
            }
            
            return clonedProblem;
        }



    template <class T_FunctionalForm>
        std::unique_ptr<dolfin::GenericLinearSolver> 
        Functional<T_FunctionalForm>::
        createSolver ()
        {
            std::string currentSolverType = parameters["current_solver_type"];
            std::string currentSolverMethod = parameters["current_solver_method"];
            std::string currentSolverPreconditioner = parameters["current_solver_preconditioner"];
            std::string desiredSolverType = parameters["desired_solver_type"];
            std::string desiredSolverMethod = parameters["desired_solver_method"];
            std::string desiredSolverPreconditioner = parameters["desired_solver_preconditioner"];
            
            if (desiredSolverType == "lu_solver")
            {
                dolfin::log (dolfin::DBG, "Creating lu_solver...");
                std::unique_ptr<dolfin::GenericLinearSolver> solver (new dolfin::LUSolver (desiredSolverMethod));
                
                dolfin::log (dolfin::DBG, "Updating parameters...");
                
                // check if solver type has changed:
                // ----- if the solver is still the same, just update the parameters
                // ----- if the parameters do not have a parameters set corresponding to the one being created but
                //       current and desired solver set are the same, that means that there is no solver at all, so
                //       add solver parameters to parameters set
                // ----- otherwise, remove current solver parameters and add new ones
                if (parameters.has_parameter_set (desiredSolverType)) 
                {
                    parameters.update (solver -> parameters);
                }
                else if (currentSolverType == desiredSolverType 
                         && currentSolverMethod == desiredSolverMethod 
                         && currentSolverPreconditioner == desiredSolverPreconditioner)
                {
                    parameters.add (solver -> parameters);
                }
                else
                {
                    parameters.remove (currentSolverType);
                    parameters.add (solver -> parameters);
                }
                
                parameters ["current_solver_type"] = desiredSolverType;
                parameters ["current_solver_method"] = desiredSolverMethod;
                parameters ["current_solver_preconditioner"] = desiredSolverPreconditioner;

                
                return solver;
            }
            else if (desiredSolverType == "krylov_solver")
            {
                dolfin::log (dolfin::DBG, "Creating krylov_solver...");
                std::unique_ptr<dolfin::GenericLinearSolver> solver (new dolfin::KrylovSolver (desiredSolverMethod, 
                                                                                               desiredSolverPreconditioner));
                
                dolfin::log (dolfin::DBG, "Updating parameters...");
                
                // check if solver type has changed:
                // ----- if the solver is still the same, just update the parameters
                // ----- if the parameters do not have a parameters set corresponding to the one being created but
                //       current and desired solver set are the same, that means that there is no solver at all, so
                //       add solver parameters to parameters set
                // ----- otherwise, remove current solver parameters and add new ones
                if (parameters.has_parameter_set (desiredSolverType)) 
                {
                    parameters.update (solver -> parameters);
                }
                else if (currentSolverType == desiredSolverType 
                         && currentSolverMethod == desiredSolverMethod 
                         && currentSolverPreconditioner == desiredSolverPreconditioner)
                {
                    parameters.add (solver -> parameters);
                }
                else
                {
                    parameters.remove (currentSolverType);
                    parameters.add (solver -> parameters);
                }
                
                parameters ["current_solver_type"] = desiredSolverType;
                parameters ["current_solver_method"] = desiredSolverMethod;
                parameters ["current_solver_preconditioner"] = desiredSolverPreconditioner;

                
                return solver;
            }
            else
            {
                dolfin::log (dolfin::DBG, "Creating solver of type \"%s\"...", desiredSolverType.c_str ());
                controlproblem::LinearSolverFactory& factory = controlproblem::LinearSolverFactory::Instance ();
                auto solver = factory.create (desiredSolverType);
                
                dolfin::log (dolfin::DBG, "Updating parameters...");
                
                // check if solver type has changed:
                // ----- if the solver is still the same, just update the parameters
                // ----- if the parameters do not have a parameters set corresponding to the one being created but
                //       current and desired solver set are the same, that means that there is no solver at all, so
                //       add solver parameters to parameters set
                // ----- otherwise, remove current solver parameters and add new ones
                if (parameters.has_parameter_set (desiredSolverType)) 
                {
                    parameters.update (solver -> parameters);
                }
                else if (currentSolverType == desiredSolverType 
                         && currentSolverMethod == desiredSolverMethod 
                         && currentSolverPreconditioner == desiredSolverPreconditioner)
                {
                    parameters.add (solver -> parameters);
                }
                else
                {
                    parameters.remove (currentSolverType);
                    parameters.add (solver -> parameters);
                }
                
                parameters ["current_solver_type"] = desiredSolverType;
                parameters ["current_solver_method"] = desiredSolverMethod;
                parameters ["current_solver_preconditioner"] = desiredSolverPreconditioner;

                
                return solver;
            }
        }
}
#endif

