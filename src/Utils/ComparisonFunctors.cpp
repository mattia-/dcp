#include <Utils/ComparisonFunctors.hpp>
#include <utility>
#include <tuple>

namespace control_problem
{
    bool less::operator() (const std::tuple <std::string, std::string, std::string>& lhs, 
                           const std::tuple <std::string, std::string, std::string>& rhs)
    {
        // if the first arguments of the two tuples are different, check which one is lower.
        // If they are the same, perform the same check on the second arguments.
        // If they still are equal, check third arguments
         
        if (std::get<0> (lhs) != std::get<0> (rhs))
        {
            return std::get<0> (lhs) < std::get<0> (rhs);
        }
        else
        {
            if (std::get<1> (lhs) != std::get<1> (rhs))
            {
                return std::get<1> (lhs) < std::get<1> (rhs);
            }
            else
            {
                return std::get<2> (lhs) < std::get<2> (rhs);
            }
        }
    }
}

