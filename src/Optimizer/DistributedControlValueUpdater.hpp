#ifndef SRC_OPTIMIZER_DISTRIBUTEDCONTROLVALUEUPDATER_HPP_INCLUDE_GUARD
#define SRC_OPTIMIZER_DISTRIBUTEDCONTROLVALUEUPDATER_HPP_INCLUDE_GUARD

#include <DifferentialProblem/CompositeDifferentialProblem.hpp>
#include <dolfin/function/GenericFunction.h>
#include <boost/shared_ptr.hpp>
#include <string>

namespace controlproblem
{
    /*! \class DistributedControlValueUpdater DistributedControlValueUpdater.hpp
     *  \brief Class to update the value of the control variable in distributed control problems.
     *  
     *  This class is a functor which can be passed to the method \c apply of any class
     *  of the \c AbstractOptimizer hierarchy, which will use it to update the value of
     *  the control parameter in the \c DifferentialProblem (also passed to the method \c apply 
     *  of the same class) as the optimization proceeds.
     */
    class DistributedControlValueUpdater
    {
        // ---------------------------------------------------------------------------------------------//

        public:
            /************************* CONSTRUCTORS ********************/
            //! Default constructor is deleted
            DistributedControlValueUpdater () = delete;
            
            //! Constructor
            /*! 
             *  Input arguments are:
             *  \param problemName string that identifies the problem (in the \c CompositeDifferentialProblem object 
             *  passed as input to \c this->operator()() ) which contains the control parameter to be updated
             *  \param coefficientType the type of the coefficient representing the control parameter inside the problem. 
             *  This will be used by the call to \c controlproblem::AbstractDifferentialProblem::setCoefficient()
             *  \param coefficientName the name of the coefficient representing the control parameter in the problem 
             *  passed as first argument
             */
            DistributedControlValueUpdater (const std::string& problemName, 
                                            const std::string& coefficientType,
                                            const std::string& coefficientName);
            
            
            /************************* OPERATORS ********************/
            //! Call operator. This will actually perform the updating of the control parameter
            /*! 
             *  Input parametes are:
             *  \param compositeProblem the problem on which to operate
             *  \param coefficientValue the new value for the control parameter identified by \c coefficientName_
             */
            void operator() (controlproblem::CompositeDifferentialProblem& compositeProblem, 
                             const boost::shared_ptr <const dolfin::GenericFunction> coefficientValue) const;
            
        // ---------------------------------------------------------------------------------------------//
        
        protected:
            //! The name identifying the problem inside the \c CompositeDifferentialProblem which contains the
            //! control parameter
            std::string problemName_;
           
            //! The type of the coefficient representing the control parameter inside the problem. This will be used
            //! by the call to \c controlproblem::AbstractDifferentialProblem::setCoefficient()
            std::string coefficientType_;
           
            //! The name of the control parameter inside the problem identified by \c problemName_
            std::string coefficientName_;
    };
}

#endif

