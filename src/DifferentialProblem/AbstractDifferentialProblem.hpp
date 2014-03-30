#ifndef SRC_DIFFERENTIALPROBLEM_ABSTRACTDIFFERENTIALPROBLEM_HPP_INCLUDE_GUARD
#define SRC_DIFFERENTIALPROBLEM_ABSTRACTDIFFERENTIALPROBLEM_HPP_INCLUDE_GUARD

#include <dolfin.h>
#include <Utils/SubdomainType.hpp>
#include <memory>

namespace control_problem
{
    /*! \class AbstractDifferentialProblem AbstractDifferentialProblem.hpp
     *  \brief Abstract base class for differential problems. 
     *         
     *  This class contains the basic interface for a differential problem to be
     *  solved with FEniCS library. It is an abstract class, it only provides the
     *  basic interface to all differential problems
     */ 
    class AbstractDifferentialProblem
    {
        // ---------------------------------------------------------------------------------------------//

        public:
            /************************* CONSTRUCTORS ********************/
            //! Default constructor
            AbstractDifferentialProblem () = delete;

            //!  Constructor with shared pointers
            /*!
             *  \param mesh the problem mesh as a const \c std::shared_ptr to \c dolfin::Mesh
             *  \param functionSpace the problem finite element space as a const \c std::shared_ptr to 
             *  \c dolfin::FunctionSpace
             *  The stored mesh's and function space's ownership will be shared between the object and the input argument.
             *  The bilinear and linear form will be created too, calling the constructor which takes the function space
             *  as input.
             */
            AbstractDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
                                         const std::shared_ptr<dolfin::FunctionSpace> functionSpace);


            //! Constructor with references
            /*!
             *  \param mesh the problem mesh as a const \c dolfin::Mesh&
             *  \param functionSpace the problem finite element space as a const \c dolfin::FunctionSpace&
             *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
             *  initialized using the \c new operator and mesh's and functionSpace's copy constructor
             *  The bilinear and linear form will be created too, calling the constructor which takes the function space
             *  as input.
             */
            AbstractDifferentialProblem (const dolfin::Mesh& mesh, 
                                         const dolfin::FunctionSpace& functionSpace);

            //! Constructor with rvalue references
            /*!
             *  \param mesh the problem mesh as a \c dolfin::Mesh&&
             *  \param functionSpace the problem finite element space as a \c dolfin::FunctionSpace&&
             *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
             *  initialized using the \c new operator and mesh's and functionSpace's move constructor
             *  The bilinear and linear form will be created too, calling the constructor which takes the function space
             *  as input.
             */
            AbstractDifferentialProblem (dolfin::Mesh&& mesh, 
                                         dolfin::FunctionSpace&& functionSpace);


            /************************* DESTRUCTOR ********************/
            //! Destructor
            /*! Default destructor, since members of the class are trivially 
             * destructible.
             * It is declared virtual so that derived classes' constructor
             * can be called on derived classes.
             * The "default-ness" is set in implementation outside of the class for compatibility with
             * gcc-4.6, which does not allow virtual members to be defaulted in class
             */
            virtual ~AbstractDifferentialProblem ();


            /********************** GETTERS ***********************/
            //! Get problem's mesh
            /*! 
             *  \return a const reference to the problem's mesh
             */
            const dolfin::Mesh& mesh () const;

            //! Get problem's finite element space
            /*! 
             *  \return a const reference to the problem's function space
             */
            const dolfin::FunctionSpace& functionSpace () const;

            //! Get const reference to the problem's i-th dirichlet boundary condition
            /*! 
             *  \param i the position in the vector storing all dirichlet BCs of the dirichlet BC object to be 
             *            returned. If called with no argument, \f$i = 0\f$ is assumed. 
             *  \return a const reference to the problem's dirichletBC in position i.
             *          No check is performed on the input value
             */
            const dolfin::DirichletBC& dirichletBC (const std::size_t& i) const;

            //! Get const reference to the problem's dirichlet boundary conditions vector
            /*! 
             *  \return a const reference to the problem's \c dirichletBC vector
             */
            const std::vector<dolfin::DirichletBC>& dirichletBCs () const;

            //! Get const reference to the problem's solution
            /*!
             *  \return a const reference to the problem's solution
             */
            const dolfin::Function& solution () const;  

            /********************** SETTERS ***********************/
            //! Set problem coefficients [1]
            /*!
             *  This method is meant to be overridden in derived classes, so that we can have a uniform interface
             *  to set coefficients in all hierarchy.
             *  Parameters are:
             *  \param coefficientType used to disambiguate between different member variables to choose which coefficient
             *  \param coefficientValue value of the coefficient
             *  \param coefficientName string identifying the coefficient
             *  to set
             */
            virtual void setCoefficient (const std::string& coefficientType, 
                                         const boost::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                         const std::string& coefficientName) = 0;

            //! Set problem coefficients [2]
            /*!
             *  This method is meant to be overridden in derived classes, so that we can have a uniform interface
             *  to set coefficients in all hierarchy.
             *  Parameters are:
             *  \param coefficientType used to disambiguate between different member variables to choose which coefficient
             *  \param coefficientValue value of the coefficient
             *  \param coefficientNumber integer identifying the coefficient
             *  to set
             */
            virtual void setCoefficient (const std::string& coefficientType,
                                         const boost::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                         const std::size_t& coefficientNumber) = 0;

            //! Set integration subdomains for the forms
            /*! 
             *  This method is meant to be overridden in derived classes, so that we can have a uniform interface
             *  to set integration subdomains in all hierarchy.
             *  Input arguments are:
             *  \param formType used to disambiguate between different member variables to choose which integration
             *  subdomain to set
             *  \param meshFunction the mesh function used to set the integration subdomains
             *  \param subdomainType the type of the subdomains, chosen among those provided by the enumeration
             *  class \c control_problem::SubdomainType
             */
            virtual void setIntegrationSubdomains (const std::string& formType,
                                                   boost::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                                   const control_problem::SubdomainType& subdomainType) = 0;

            //! Add Dirichlet boundary condition to the problem [1]
            /*!
             *  \param dirichletCondition a const reference to the dirichlet boundary condition to be added to the problem
             */
            virtual void addDirichletBC (const dolfin::DirichletBC& dirichletCondition);

            //! Add Dirichlet boundary condition to the problem [2]
            /*!
             *  \param dirichletCondition a rvalue reference to the dirichlet boundary condition to be added to the problem
             */
            virtual void addDirichletBC (dolfin::DirichletBC&& dirichletCondition);

            //! Remove Dirichlet boundary condition with given position
            /*!
             *  \param i the position in the vector of the boundary condition to be removed.
             *            If i is greater than the size of the vector, nothing is removed.
             */
            virtual void removeDirichletBC (const std::vector<dolfin::DirichletBC>::iterator& i);
            
            //! This method is meant to be overridden only if needed in the derived classes. It checks for possible
            //! private members to update (e.g. the solver, if the solver method string in the parameters
            //! was changed but the solver itself was not) and performs the updating. In this class, it is implemented
            //! as an empty funciton, so there is no need to override it if the derived class has no need for this
            //! method
            virtual void update ();

            /********************** METHODS ***********************/

            //! Solve method
            /*!
             * Solves differential problem storing the solution in the private member \c solution_.
             * It is a pure virtual method that needs to be overridden
             * in any concrete instance of the class
             */
            virtual void solve () = 0;
            
            //! Clone method
            /*!
             *  \return a pointer to a \c control_problem::AbstractDifferentialProblem containing a copy of the object on 
             *  which it is called. 
             */
            virtual control_problem::AbstractDifferentialProblem* clone () const = 0;

            /********************** VARIABLES ***********************/
            //! the problem parameters
            dolfin::Parameters parameters;
            
            // ---------------------------------------------------------------------------------------------//

        protected:
            //! The problem mesh
            /*! 
             *  Stored as a \c shared_ptr because it may be common to more than 
             *  one problem
             */
            std::shared_ptr<dolfin::Mesh> mesh_;

            //! The problem finite element space
            /*! 
             *  Stored as a \c shared_ptr because it may be common to more than 
             *  one problem
             */
            std::shared_ptr<dolfin::FunctionSpace> functionSpace_;

            //! The Dirichlet's boundary conditions vector
            std::vector<dolfin::DirichletBC> dirichletBCs_;

            //! Solution of the differential problem
            dolfin::Function solution_;

            // ---------------------------------------------------------------------------------------------//

        private:

    };

}
#endif
