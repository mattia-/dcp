#ifndef SRC_DIFFERENTIALPROBLEM_LINEARDIFFERENTIALPROBLEM_HPP_INCLUDE_GUARD
#define SRC_DIFFERENTIALPROBLEM_LINEARDIFFERENTIALPROBLEM_HPP_INCLUDE_GUARD

#include <dolfin.h>
#include <vector>
#include <string>
#include <memory>
#include <DifferentialProblem/AbstractDifferentialProblem.hpp>
#include <Factory/LinearSolverFactory.hpp>
#include <Utils/SubdomainType.hpp>

namespace controlproblem
{
    /*! \class LinearDifferentialProblem LinearDifferentialProblem.hpp
     *  \brief Class for linear differential problems.
     *
     *  This class represents problem of the form
     *  \f[
     *      \mbox{Find } u \in V : a \left(u, v\right) = F \left(v\right) \ \forall\,v\,\in\,V
     *  \f]
     *  with \f$ a \left(u, v\right) : V \times V \rightarrow \mathds{R}\f$ bilinear form on \f$V\f$
     *  and \f$ L \left(v\right) : V \rightarrow \mathds{R} \f$ linear form on the same space.
     *  
     *  It inherits publicly from \c AbstractDifferentialProblem
     *  and it extends its functionalities to a concrete differential
     *  problem.
     *  Template arguments are:
     *  \arg T_BilinearForm the bilinear form type
     *  \arg T_LinearForm the linear form type
     *  \arg T_LinearSolverFactory the type of the factory that creates the linear solver. By default, it is set
     *  to \c controlproblem::LinearSolverFactory
     */

    template <class T_BilinearForm_, class T_LinearForm_, class T_LinearSolverFactory_ = controlproblem::LinearSolverFactory>
        class LinearDifferentialProblem : public AbstractDifferentialProblem
        {
            // ---------------------------------------------------------------------------------------------//  

            public:
                /******************* TYPEDEFS *******************/
                typedef T_BilinearForm_        T_BilinearForm;
                typedef T_LinearForm_          T_LinearForm;
                typedef T_LinearSolverFactory_ T_LinearSolverFactory;
                
                
                /******************* CONSTRUCTORS *******************/
                //! Default constructor is deleted. The class is not default constructable.
                LinearDifferentialProblem () = delete;

                //!  Constructor with shared pointers [1]
                /*!
                 *  \param mesh the problem mesh as a \c const \c std::shared_ptr to \c dolfin::Mesh
                 *  \param functionSpace the problem finite element space as a \c const \c std::shared_ptr to 
                 *  \c dolfin::FunctionSpace
                 *  \param solverType the type of the solver. Default: \c lu_solver
                 *  \param solverMethod the method of the solver. Possible values depend on the solver type (see dolfin
                 *  documentation, or use method \c list_<solverType>_methods). Default value: \c default
                 *  \param solverPreconditioner the preconditioner to be used. It is not used for \c lu_solvers. Default
                 *  value: \c default
                 *  The stored mesh's and function space's ownership will be shared between the object and the input argument.
                 *  The bilinear and linear form will be created too, calling the constructor which takes the function space
                 *  as input.
                 */
                LinearDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
                                           const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
                                           const std::string& solverType = "lu_solver",
                                           const std::string& solverMethod = "default",
                                           const std::string& solverPreconditioner = "default");
                

                //! Constructor with references [1]
                /*!
                 *  \param mesh the problem mesh as a \c const \c dolfin::Mesh&
                 *  \param functionSpace the problem finite element space as a \c const \c dolfin::FunctionSpace&
                 *  \param solverType the type of the solver. Default: \c lu_solver
                 *  \param solverMethod the method of the solver. Possible values depend on the solver type (see dolfin
                 *  documentation, or use method \c list_<solverType>_methods). Default value: \c default
                 *  \param solverPreconditioner the preconditioner to be used. It is not used for lu solvers. Default
                 *  value: \c default
                 *  The stored mesh's and function space's ownership will be unique to the object, since the protected
                 *  member variables are 
                 *  initialized using the \c new operator and mesh's and functionSpace's copy constructor.
                 *  The bilinear and linear form will be created too, calling the constructor which takes the function space
                 *  as input.
                 */
                LinearDifferentialProblem (const dolfin::Mesh& mesh, 
                                           const dolfin::FunctionSpace& functionSpace,
                                           const std::string& solverType = "lu_solver",
                                           const std::string& solverMethod = "default",
                                           const std::string& solverPreconditioner = "default");

                //! Constructor with rvalue references [1]
                /*!
                 *  \param mesh the problem mesh as a \c dolfin::Mesh&&
                 *  \param functionSpace the problem finite element space as a \c dolfin::FunctionSpace&&
                 *  \param solverType the type of the solver. Default: \c lu_solver
                 *  \param solverMethod the method of the solver. Possible values depend on the solver type (see dolfin
                 *  documentation, or use method \c list_<solverType>_methods). Default value: \c default
                 *  \param solverPreconditioner the preconditioner to be used. It is not used for lu solvers. Default
                 *  value: \c default
                 *  The stored mesh's and function space's ownership will be unique to the object, since the protected
                 *  member variables are 
                 *  initialized using the \c new operator and mesh's and functionSpace's move constructor
                 *  The bilinear and linear form will be created too, calling the constructor which takes the function space
                 *  as input.
                 */
                LinearDifferentialProblem (dolfin::Mesh&& mesh, 
                                           dolfin::FunctionSpace&& functionSpace,
                                           const std::string& solverType = "lu_solver",
                                           const std::string& solverMethod = "default",
                                           const std::string& solverPreconditioner = "default");

                
                //!  Constructor with shared pointers [2]
                /*!
                 *  \param mesh the problem mesh as a \c const \c std::shared_ptr to \c dolfin::Mesh
                 *  \param functionSpace the problem finite element space as a \c const \c std::shared_ptr to 
                 *  \c dolfin::FunctionSpace
                 *  \param bilinearForm a \c const reference to the problem's bilinear form
                 *  \param linearForm a \c const reference to the problem's linear form
                 *  \param solverType the type of the solver. Default: \c lu_solver
                 *  \param solverMethod the method of the solver. Possible values depend on the solver type (see dolfin
                 *  documentation, or use method \c list_<solverType>_methods. Default value: \c default
                 *  \param solverPreconditioner the preconditioner to be used. It is not used for lu solvers. Default
                 *  value: \c default
                 *  The stored mesh's and function space's ownership will be shared between the object and the input argument.
                 *  The bilinear and linear form will be created too, calling the constructor which takes the function space
                 *  as input.
                 */
                LinearDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
                                           const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
                                           const T_BilinearForm& bilinearForm,
                                           const T_LinearForm& linearForm,
                                           const std::string& solverType = "lu_solver",
                                           const std::string& solverMethod = "default",
                                           const std::string& solverPreconditioner = "default");

                //! Constructor with references [2]
                /*!
                 *  \param mesh the problem mesh as a \c const \c dolfin::Mesh&
                 *  \param functionSpace the problem finite element space as a \c const \c dolfin::FunctionSpace&
                 *  \param bilinearForm a \c const reference to the problem's bilinear form
                 *  \param linearForm a \c const reference to the problem's linear form
                 *  \param solverType the type of the solver. Default: \c lu_solver
                 *  \param solverMethod the method of the solver. Possible values depend on the solver type (see dolfin
                 *  documentation, or use method \c list_<solverType>_methods. Default value: \c default
                 *  \param solverPreconditioner the preconditioner to be used. It is not used for lu solvers. Default
                 *  value: \c default
                 *  The stored mesh's and function space's ownership will be unique to the object, since the protected
                 *  member variables are 
                 *  initialized using the \c new operator and mesh's and functionSpace's copy constructor
                 *  The bilinear and linear form will be created too, calling the constructor which takes the function space
                 *  as input.
                 */
                LinearDifferentialProblem (const dolfin::Mesh& mesh, 
                                           const dolfin::FunctionSpace& functionSpace,
                                           const T_BilinearForm& bilinearForm,
                                           const T_LinearForm& linearForm,
                                           const std::string& solverType = "lu_solver",
                                           const std::string& solverMethod = "default",
                                           const std::string& solverPreconditioner = "default");

                //! Constructor with rvalue references [3]
                /*!
                 *  \param mesh the problem mesh as a \c dolfin::Mesh&&
                 *  \param functionSpace the problem finite element space as a \c dolfin::FunctionSpace&&
                 *  \param bilinearForm a rvalue reference to the problem's bilinear form
                 *  \param linearForm a rvalue reference to the problem's linear form
                 *  \param solverType the type of the solver. Default: \c lu_solver
                 *  \param solverMethod the method of the solver. Possible values depend on the solver type (see dolfin
                 *  documentation, or use method \c list_<solverType>_methods. Default value: \c default
                 *  \param solverPreconditioner the preconditioner to be used. It is not used for lu solvers. Default
                 *  value: \c default
                 *  The stored mesh's and function space's ownership will be unique to the object, since the protected
                 *  member variables are 
                 *  initialized using the \c new operator and mesh's and functionSpace's move constructor
                 *  The bilinear and linear form will be created too, calling the constructor which takes the function space
                 *  as input.
                 */
                LinearDifferentialProblem (dolfin::Mesh&& mesh, 
                                           dolfin::FunctionSpace&& functionSpace,
                                           T_BilinearForm&& bilinearForm,
                                           T_LinearForm&& linearForm,
                                           const std::string& solverType = "lu_solver",
                                           const std::string& solverMethod = "default",
                                           const std::string& solverPreconditioner = "default");

                /******************* DESTRUCTOR *******************/
                
                //! Destructor
                /*! 
                 *  Default destructor, since members of the class are trivially 
                 *  destructible.
                 *  It is declared virtual so that derived classes' constructor
                 *  can be called on derived classes.
                 *  The "default-ness" is set in implementation outside of the class for compatibility with
                 *  \c gcc-4.6, which does not allow virtual members to be defaulted in class
                 */
                virtual ~LinearDifferentialProblem ();

                
                /******************* GETTERS *******************/
                //! Get const reference to the problem's linear form
                /*! 
                 *  \return a const reference to the problem's linear form
                 */
                const T_BilinearForm& bilinearForm () const;

                //! Get const reference to the problem's linear form
                /*! 
                 *  \return a const reference to the problem's linear form
                 */
                const T_LinearForm& linearForm () const;

                //! Get const reference to the problem's linear operator
                /*!
                 *  \return a const reference to the problem's linear operator, which
                 *  is a \c dolfin::Matrix
                 */
                const dolfin::Matrix& linearOperator () const;

                //! Get const reference to the problem's right hand side
                /*!
                 *  \return a const reference to the problem's right hand side, which
                 *  is a \c dolfin::Vector
                 */
                const dolfin::Vector& rhs () const;


                /******************* SETTERS *******************/
                
                //! Set coefficient [1]. Override of virtual function in \c AbstractDifferentialProblem.
                /*!
                 *  Possible values for \c coefficientType are:
                 *  \li \c bilinear_form to set the coefficient in the bilinear form
                 *  \li \c linear_form to set the coefficient in the linear form
                 *  
                 *  See \c AbstractDifferentialProblem documentation for more details on the function
                 */
                virtual void setCoefficient (const std::string& coefficientType, 
                                             const boost::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                             const std::string& coefficientName);

                //! Set coefficient [2]. Override of virtual function in \c AbstractDifferentialProblem.
                /*!
                 *  Possible values for \c coefficientType are:
                 *  \li \c bilinear_form to set the coefficient in the bilinear form
                 *  \li \c linear_form to set the coefficient in the linear form
                 *  
                 *  See \c AbstractDifferentialProblem documentation for more details on the function
                 */
                virtual void setCoefficient (const std::string& coefficientType,
                                             const boost::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                             const std::size_t& coefficientNumber);

                //! Set integration subdomains for the forms. Override of virtual function in \c AbstractDifferentialProblem
                /*! 
                 *  Possible values for \c coefficientType are:
                 *  \li \c bilinear_form to set the integration subdomain in the bilinear form
                 *  \li \c linear_form to set the integration subdomain in the linear form
                 *  
                 *  See \c AbstractDifferentialProblem documentation for more details on the function
                 */
                virtual void setIntegrationSubdomains (const std::string& formType,
                                                       boost::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                                       const controlproblem::SubdomainType& subdomainType);

                //! Add Dirichlet boundary condition to the problem [1]. Overrides method in \c AbstractDifferentialProblem
                /*!
                 *  This method adds to the base class method the setting of parameter \c system_is_assembled to \c false.
                 *  \param dirichletCondition a const reference to the dirichlet boundary condition to be added to the problem
                 */
                virtual void addDirichletBC (const dolfin::DirichletBC& dirichletCondition);

                //! Add Dirichlet boundary condition to the problem [2]. Overrides method in \c AbstractDifferentialProblem
                /*!
                 *  This method adds to the base class method the setting of parameter \c system_is_assembled to \c false.
                 *  \param dirichletCondition a rvalue reference to the dirichlet boundary condition to be added to the problem
                 */
                virtual void addDirichletBC (dolfin::DirichletBC&& dirichletCondition);

                //! Remove Dirichlet boundary condition with given position. Overrides method in \c AbstractDifferentialProblem
                /*!
                 *  This method adds to the base class method the setting of parameter \c system_is_assembled to \c false.
                 *  \param i the position in the vector of the boundary condition to be removed.
                 *            If i is greater than the size of the vector, nothing is removed.
                 */
                virtual void removeDirichletBC (const std::vector<dolfin::DirichletBC>::iterator& i);
                
                //! Method to update class members. It checks for differences between desired and current solver parameters
                //! and creates a new solver, setting also the proper parameters
                virtual void update ();

                
                /******************* METHODS *******************/

                //! Solve problem
                /*!
                 *  This method solves the problem defined. It uses the private members' value to set the problem and then
                 *  stores the solution in the private member \c solution_. Note that it checks the problem's member
                 *  \c parameters to decided whether problem's matrix and vector should be reassembled and if
                 *  the values of the parameters \c desired_solver_type, \c desired_solver_method and 
                 *  \c desired_solver_preconditioner match the values of \c current_solver_type, \c current_solver_method
                 *  and \c current_solver_preconditioner. If they differ, it calls \c createSolver()
                 */
                virtual void solve ();

                //! Solve problem specifying flag
                /*!
                 *  \param mustReassemble set it to \c true if the system operators (matrix and right hand side vector)
                 *         should be reassembled. It is \c false by default.
                 */
                void solve (const bool& mustReassemble);

                //! Clone method. Overrides method in \c AbstractDifferentialProblem
                /*!
                 *  Note that it uses variable \c clone_method in \c parameters to decide which kind of cloning to
                 *  perform. Values for such variable can be either \c deep_clone or \c shallow_clone. The first means
                 *  that the new object is created calling the constructor that takes a mesh and a function space as 
                 *  input, thus creating a copy of such objects and returning a completely independent cloned object. 
                 *  The second cloning type calls the constructor that takes shared pointers as input: the mesh and
                 *  the function space are not copied but shared between the current object and its clone. 
                 *  The default value for \clone_method is \shallow_clone
                 */
                virtual controlproblem::LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>*
                    clone () const;

                // ---------------------------------------------------------------------------------------------//

            protected:
                //! Creates a linear solver object of type passed as input. 
                /*!
                 *  The solver will be created using the parameters \c desired_solver_type, \c desired_solver_method
                 *  and \c desired_solver_preconditioner set in the protected member \c parameters. It will
                 *  also set the parameters \c current_solver_type, \c current_solver_method and 
                 *  \c current_solver_preconditioner in the same set of parameters substituting the current values with 
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
                 *  We use a \c boost::shared_ptr for compatibility with FEniCS, version 1.3.0.
                 *  Note that there is no write access to this variable from outside the class,
                 *  so it is guaranteed to be a unique shared pointer
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

    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        LinearDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
                                   const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
                                   const std::string& solverType,
                                   const std::string& solverMethod,
                                   const std::string& solverPreconditioner) :
            AbstractDifferentialProblem (mesh, functionSpace),
            bilinearForm_ (*functionSpace_, *functionSpace_),
            linearForm_ (*functionSpace_),
            solver_ (nullptr),
            problemMatrix_ (new dolfin::Matrix),
            rhsVector_ ()
        { 
            dolfin::begin (dolfin::DBG, "Building LinearDifferentialProblem...");
            
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
            parameters.add ("clone_method", "shallow_clone");
            
            dolfin::begin (dolfin::DBG, "Creating solver...");
            solver_ = createSolver ();
            dolfin::end ();
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "LinearDifferentialProblem created");
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        LinearDifferentialProblem (const dolfin::Mesh& mesh, 
                                   const dolfin::FunctionSpace& functionSpace,
                                   const std::string& solverType,
                                   const std::string& solverMethod,
                                   const std::string& solverPreconditioner) :
            AbstractDifferentialProblem (mesh, functionSpace),
            bilinearForm_ (*functionSpace_, *functionSpace_),
            linearForm_ (*functionSpace_),
            solver_ (nullptr),
            problemMatrix_ (new dolfin::Matrix),
            rhsVector_ ()
        { 
            dolfin::begin (dolfin::DBG, "Building LinearDifferentialProblem...");
            
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
            parameters.add ("clone_method", "shallow_clone");
            
            dolfin::begin (dolfin::DBG, "Creating solver...");
            solver_ = createSolver ();
            dolfin::end ();
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "LinearDifferentialProblem created");
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        LinearDifferentialProblem (dolfin::Mesh&& mesh, 
                                   dolfin::FunctionSpace&& functionSpace,
                                   const std::string& solverType,
                                   const std::string& solverMethod,
                                   const std::string& solverPreconditioner) :
            AbstractDifferentialProblem (mesh, functionSpace),
            bilinearForm_ (*functionSpace_, *functionSpace_),
            linearForm_ (*functionSpace_),
            solver_ (nullptr),
            problemMatrix_ (new dolfin::Matrix),
            rhsVector_ ()
        { 
            dolfin::begin (dolfin::DBG, "Building LinearDifferentialProblem...");
            
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
            parameters.add ("clone_method", "shallow_clone");
        
            dolfin::begin (dolfin::DBG, "Creating solver...");
            solver_ = createSolver ();
            dolfin::end ();
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "LinearDifferentialProblem created");
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        LinearDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
                                   const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
                                   const T_BilinearForm& bilinearForm,
                                   const T_LinearForm& linearForm,
                                   const std::string& solverType,
                                   const std::string& solverMethod,
                                   const std::string& solverPreconditioner) :
            AbstractDifferentialProblem (mesh, functionSpace),
            bilinearForm_ (bilinearForm),
            linearForm_ (linearForm),
            solver_ (nullptr),
            problemMatrix_ (new dolfin::Matrix),
            rhsVector_ ()
        { 
            dolfin::begin (dolfin::DBG, "Building LinearDifferentialProblem...");
            
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
            parameters.add ("clone_method", "shallow_clone");
            
            dolfin::begin (dolfin::DBG, "Creating solver...");
            solver_ = createSolver ();
            dolfin::end ();
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "LinearDifferentialProblem created");
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        LinearDifferentialProblem (const dolfin::Mesh& mesh, 
                                   const dolfin::FunctionSpace& functionSpace,
                                   const T_BilinearForm& bilinearForm,
                                   const T_LinearForm& linearForm,
                                   const std::string& solverType,
                                   const std::string& solverMethod,
                                   const std::string& solverPreconditioner) :
            AbstractDifferentialProblem (mesh, functionSpace),
            bilinearForm_ (bilinearForm),
            linearForm_ (linearForm),
            solver_ (nullptr),
            problemMatrix_ (new dolfin::Matrix),
            rhsVector_ ()
        { 
            dolfin::begin (dolfin::DBG, "Building LinearDifferentialProblem...");
            
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
            parameters.add ("clone_method", "shallow_clone");
            
            dolfin::begin (dolfin::DBG, "Creating solver...");
            solver_ = createSolver ();
            dolfin::end ();
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "LinearDifferentialProblem created");
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        LinearDifferentialProblem (dolfin::Mesh&& mesh, 
                                   dolfin::FunctionSpace&& functionSpace,
                                   T_BilinearForm&& bilinearForm,
                                   T_LinearForm&& linearForm,
                                   const std::string& solverType,
                                   const std::string& solverMethod,
                                   const std::string& solverPreconditioner) :
            AbstractDifferentialProblem (mesh, functionSpace),
            bilinearForm_ (std::move (bilinearForm)),
            linearForm_ (std::move (linearForm)),
            solver_ (nullptr),
            problemMatrix_ (new dolfin::Matrix),
            rhsVector_ ()
        {
            dolfin::begin (dolfin::DBG, "Building LinearDifferentialProblem...");
            
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
            parameters.add ("clone_method", "shallow_clone");
            
            dolfin::begin (dolfin::DBG, "Creating solver...");
            solver_ = createSolver ();
            dolfin::end ();
            
            dolfin::end ();
            
            dolfin::log (dolfin::DBG, "LinearDifferentialProblem created");
        }

    
    /***************** DESTRUCTOR ******************/

    // this is done for compatibility with gcc-4.6, which doesn't allow virtual members to be defualted in class body
    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        ~LinearDifferentialProblem () = default;
    
    


    /******************* GETTERS *******************/

    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        const T_BilinearForm& LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        bilinearForm () const
        {
            return bilinearForm_;
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        const T_LinearForm& LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        linearForm () const
        {
            return linearForm_;
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        const dolfin::Matrix& LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        linearOperator () const
        {
            return *problemMatrix_;
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        const dolfin::Vector& LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        rhs () const
        {
            return rhsVector_;
        }



    /******************* SETTERS *******************/

    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
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

    

    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
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



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
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



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        addDirichletBC (const dolfin::DirichletBC& dirichletCondition)
        {
            dolfin::log (dolfin::DBG, "Adding dirichlet boundary condition to boundary conditions vector...");
            dirichletBCs_.emplace_back (dirichletCondition);
            parameters ["system_is_assembled"] = false;
        }

    

    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        addDirichletBC (dolfin::DirichletBC&& dirichletCondition)
        {
            dolfin::log (dolfin::DBG, "Adding dirichlet boundary condition to boundary conditions vector...");
            dirichletBCs_.emplace_back (dirichletCondition);
            parameters ["system_is_assembled"] = false;
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        removeDirichletBC (const std::vector<dolfin::DirichletBC>::iterator& i)
        {
            dolfin::log (dolfin::DBG, "Removing dirichlet boundary condition from boundary conditions vector...");
            dirichletBCs_.erase (i);
            parameters ["system_is_assembled"] = false;
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        update ()
        {
            // define auxiliary string variables
            std::string desiredSolverType = parameters ["desired_solver_type"];
            std::string desiredSolverMethod = parameters ["desired_solver_method"];
            std::string desiredSolverPreconditioner = parameters ["desired_solver_preconditioner"];
            
            std::string currentSolverType = parameters ["current_solver_type"];
            std::string currentSolverMethod = parameters ["current_solver_method"];
            std::string currentSolverPreconditioner = parameters ["current_solver_preconditioner"];
            
            bool needsSolverUpdating = desiredSolverType != currentSolverType 
                                       || desiredSolverMethod != currentSolverMethod 
                                       || desiredSolverPreconditioner != currentSolverPreconditioner;
            
            if (needsSolverUpdating)
            {
                dolfin::begin (dolfin::DBG, "Updating solver...");
                solver_ = createSolver ();
                dolfin::end ();
                
                parameters ["system_is_assembled"] = false;
            }
        }



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        solve () 
        {
            update ();
            
            // define auxiliary string variables
            bool systemIsAssembled = parameters ["system_is_assembled"];
            bool forceReassembleSystem = parameters ["force_reassemble_system"];
            bool needsReassembling = !systemIsAssembled || forceReassembleSystem;
            
            if (needsReassembling)
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



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
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
    
    
    
    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        controlproblem::LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>*
        LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
        clone () const
        {
            dolfin::begin (dolfin::DBG, "Cloning object...");
            
            if (dolfin::get_log_level () > dolfin::DBG)
            {
                dolfin::end ();
            }
            
            std::string cloneMethod = parameters["clone_method"];
            
            dolfin::log (dolfin::DBG, "Clone method: %s", cloneMethod.c_str ());
            dolfin::log (dolfin::DBG, "Creating new object of type LinearDifferentialProblem...");
            
            // create new object
            controlproblem::LinearDifferentialProblem <T_BilinearForm, T_LinearForm, T_LinearSolverFactory>* clonedProblem;
            if (cloneMethod == "shallow_clone")
            {
                clonedProblem = 
                    new controlproblem::LinearDifferentialProblem <T_BilinearForm, T_LinearForm, T_LinearSolverFactory> 
                    (this->mesh_,
                     this->functionSpace_,
                     this->bilinearForm_, 
                     this->linearForm_
                    );
            }
            else if (cloneMethod == "deep_clone")
            {
                clonedProblem =
                    new controlproblem::LinearDifferentialProblem <T_BilinearForm, T_LinearForm, T_LinearSolverFactory> 
                    (*(this->mesh_),
                     *(this->functionSpace_),
                       this->bilinearForm_, 
                       this->linearForm_
                    );
            }
            else
            {
                dolfin::error ("Cannot clone linear differential problem. Unknown clone method: \"%s\"",
                               cloneMethod.c_str ());
            }
            
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



    template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
        std::unique_ptr<dolfin::GenericLinearSolver> 
        LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
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
