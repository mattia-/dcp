#include <Optimizer/DistributedControlValueUpdater.hpp>
#include <DifferentialProblem/AbstractDifferentialProblem.hpp>

namespace DCP
{
    /************************* CONSTRUCTORS ********************/
    DistributedControlValueUpdater::DistributedControlValueUpdater (const std::string& problemName, 
                                                                    const std::string& coefficientType,
                                                                    const std::string& coefficientName) :
        problemName_ (problemName),
        coefficientType_ (coefficientType),
        coefficientName_ (coefficientName)
    {   }
     


    /************************* OPERATORS ********************/
    void DistributedControlValueUpdater::operator() (DCP::CompositeDifferentialProblem& compositeProblem, 
                                                     const boost::shared_ptr <const dolfin::GenericFunction> coefficientValue) const
    {
        DCP::AbstractDifferentialProblem& problem = compositeProblem [problemName_];
        
        problem.setCoefficient (coefficientType_, coefficientValue, coefficientName_);
    }
}
