#include <Optimizer/NeumannControlValueUpdater.hpp>
#include <DifferentialProblem/AbstractDifferentialProblem.hpp>

namespace controlproblem
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
    void NeumannControlValueUpdater::operator() (controlproblem::CompositeDifferentialProblem& compositeProblem, 
                                                 const boost::shared_ptr <const dolfin::GenericFunction> coefficientValue) const
{
        controlproblem::AbstractDifferentialProblem& problem = compositeProblem [problemName_];

        problem.setCoefficient (coefficientType_, coefficientValue, coefficientName_);
    }
}
