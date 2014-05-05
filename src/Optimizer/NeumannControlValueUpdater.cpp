#include <Optimizer/NeumannControlValueUpdater.hpp>
#include <DifferentialProblem/AbstractDifferentialProblem.hpp>

namespace DCP
{
    /************************* CONSTRUCTORS ********************/
    NeumannControlValueUpdater::NeumannControlValueUpdater (const std::string& problemName, 
                                                            const std::string& coefficientType,
                                                            const std::string& coefficientName) :
        problemName_ (problemName),
        coefficientType_ (coefficientType),
        coefficientName_ (coefficientName)
    {   }



    /************************* OPERATORS ********************/
    void NeumannControlValueUpdater::operator() (DCP::CompositeDifferentialProblem& compositeProblem, 
                                                 const boost::shared_ptr <const dolfin::GenericFunction> coefficientValue) const
    {
        DCP::AbstractDifferentialProblem& problem = compositeProblem [problemName_];

        problem.setCoefficient (coefficientType_, coefficientValue, coefficientName_);
    }
}
