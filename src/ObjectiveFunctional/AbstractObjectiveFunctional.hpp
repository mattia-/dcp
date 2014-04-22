#ifndef SRC_OBJECTIVEFUNCTIONAL_ABSTRACTOBJECTIVEFUNCTIONAL_HPP_INCLUDE_GUARD
#define SRC_OBJECTIVEFUNCTIONAL_ABSTRACTOBJECTIVEFUNCTIONAL_HPP_INCLUDE_GUARD

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/common/Array.h>
#include <dolfin/fem/Form.h>
#include <boost/shared_ptr.hpp>
#include <string>
#include <ObjectiveFunctional/VariableExpression.hpp>
#include <DifferentialProblem/SubdomainType.hpp>

namespace controlproblem
{
    /*! \class AbstractObjectiveFunctional AbstractObjectiveFunctional.hpp
     *  \brief Abstract base class for objective functionals.
     *
     *  This class contains the basic interface for a generic objective functional. The class is pure virtual 
     *  and it is intended to be use in order to apply polymorphism for concrete instances of the derived 
     *  \c ObjectiveFunctional class, which is derived for this one
     */

    class AbstractObjectiveFunctional
    {
        // ---------------------------------------------------------------------------------------------//  

        public:
            /******************* CONSTRUCTORS *******************/
            //! Default constructor is deleted
            AbstractObjectiveFunctional () = delete;

            //! Constructor with \c shared_ptr
            AbstractObjectiveFunctional (const boost::shared_ptr <const dolfin::Mesh> mesh);
            
            //! Constructor with <tt> const reference </tt>
            AbstractObjectiveFunctional (const dolfin::Mesh& mesh);


            /******************* DESTRUCTOR *******************/

            //! Default destructor
            virtual ~AbstractObjectiveFunctional () {};

            
            /******************* GETTERS *******************/
            //! Get const reference to the mesh
            /*! 
             *  \return a const reference to the mesh 
             */
            virtual const dolfin::Mesh& mesh () const;

            //! Get const reference to the functional
            /*! 
             *  \return a const reference to the functional form
             */
            virtual const dolfin::Form& functional () const = 0;

            //! Get const reference to the functional gradient
            /*! 
             *  \return a const reference to the functional gradient
             */
            virtual const dolfin::Expression& gradient () const = 0;


            /******************* SETTERS *******************/
            //! Set coefficients for the protected member variables
            /*!
             *  This method is meant to be overridden in derived classes.
             *  Parameters are:
             *  \param coefficientType used to disambiguate between different member variables to choose which coefficient
             *  to set
             *  \param coefficientValue value of the coefficient
             *  \param coefficientName string identifying the coefficient to set
             */
            virtual void setCoefficient (const std::string& coefficientType, 
                                         const boost::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                         const std::string& coefficientName) = 0;

            //! Set integration subdomains for the protected member variable that represent the functional (which must
            //! be declared in the derived class)
            /*! 
             *  This method is meant to be overridden in derived classes.
             *  Input arguments are:
             *  \param meshFunction the mesh function used to set the integration subdomains
             *  \param subdomainType the type of the subdomains, chosen among those provided by the enumeration
             *  class \c controlproblem::SubdomainType
             */
            virtual void setIntegrationSubdomains (boost::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                                   const controlproblem::SubdomainType& subdomainType) = 0;

            /******************* METHODS *******************/
            //! Evaluate the stored functional
            /*!
             *  \return a double containing the evaluation of the functional on the mesh
             */
            virtual double evaluateFunctional () const = 0;

            //! Evaluate the stored functional gradient at given point
            /*!
             *  Input arguments are:
             *  \param values the values at the point
             *  \param x the coordinates of the point
             *  \param cell the cell which contains the given point
             *  
             *  \return a double containing the evaluation of the functional on the mesh
             */
            virtual void evaluateGradient (dolfin::Array<double>& values, 
                                           const dolfin::Array<double>& x, 
                                           const ufc::cell& cell) const = 0;

            //! Evaluate at given point in given cell
            /*!
             *  Input arguments are:
             *  \param values the values at the point
             *  \param x the coordinates of the point
             *  
             *  \return a double containing the evaluation of the functional on the mesh
             */
            virtual void evaluateGradient (dolfin::Array<double>& values, 
                                           const dolfin::Array<double>& x) const = 0;

            
            // ---------------------------------------------------------------------------------------------//

        protected:
            //! The mesh over which the functional is defined
            boost::shared_ptr <const dolfin::Mesh> mesh_;
            
            
            // ---------------------------------------------------------------------------------------------//

        private:
    };
}
#endif
