#include <ObjectiveFunctional/VariableExpression.hpp>
#include <utility>
#include <memory>

namespace DCP
{
    /******************* CONSTRUCTORS *******************/
    VariableExpression::VariableExpression (std::size_t dim) : 
        Expression (dim),
        variables_ ()
    { 
        dolfin::log (dolfin::DBG, "VariableExpression object created");
    }
    



    VariableExpression::VariableExpression (std::size_t dim0, std::size_t dim1) : 
        Expression (dim0, dim1),
        variables_ ()
    { 
        dolfin::log (dolfin::DBG, "VariableExpression object created");
    }
    


    VariableExpression::VariableExpression (std::vector<std::size_t> value_shape) : 
        Expression (value_shape),
        variables_ ()
    { 
        dolfin::log (dolfin::DBG, "VariableExpression object created");
    }

    

    VariableExpression::
    VariableExpression (const std::map <std::string, boost::shared_ptr <const dolfin::GenericFunction>>& variables) : 
        Expression (),
        variables_ (variables)
    {
        dolfin::log (dolfin::DBG, "VariableExpression object created");
    }
        
           

    VariableExpression::
    VariableExpression (std::size_t dim,
                        const std::map <std::string, boost::shared_ptr <const dolfin::GenericFunction>>& variables) :
        Expression (dim), 
        variables_ (variables)
    {
        dolfin::log (dolfin::DBG, "VariableExpression object created");
    }

    

    VariableExpression::
    VariableExpression (std::size_t dim0, 
                        std::size_t dim1,
                        const std::map <std::string, boost::shared_ptr <const dolfin::GenericFunction>>& variables) : 
        Expression (dim0, dim1),
        variables_ (variables)
    {
        dolfin::log (dolfin::DBG, "VariableExpression object created");
    }

    

    VariableExpression::
    VariableExpression (std::vector<std::size_t> value_shape,
                        const std::map <std::string, boost::shared_ptr <const dolfin::GenericFunction>>& variables) : 
        Expression (value_shape),
        variables_ (variables)
    {
        dolfin::log (dolfin::DBG, "VariableExpression object created");
    }

    

    /******************* SETTERS *******************/
    void VariableExpression::setCoefficient (const std::string& variableName, 
                                             const boost::shared_ptr <const dolfin::GenericFunction> value)
    {
        dolfin::log (dolfin::DBG, "Inserting variable in map...");
        
        variables_ [variableName] = value;
    }


    
    /******************* GETTERS *******************/
    const dolfin::Function& VariableExpression::function (const std::string& variableName) const
    {
        auto variable = variables_.find (variableName);
        return *(boost::dynamic_pointer_cast<const dolfin::Function> (variable -> second));
    }

    

    const dolfin::Expression& VariableExpression::expression (const std::string& variableName) const
    {
        auto variable = variables_.find (variableName);
        return *(boost::dynamic_pointer_cast<const dolfin::Expression> (variable -> second));
    }
        
    

    /******************* METHODS *******************/
    void VariableExpression::eval (dolfin::Array<double>& values, const dolfin::Array<double>& x, const ufc::cell& cell) const
    {
        // Redirect to simple eval
        eval(values, x); 
    }



    void VariableExpression::eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    {
        dolfin::dolfin_error("VariableExpression.cpp",
                             "evaluate expression",
                             "Missing eval() function (must be overloaded)");
    }



    void VariableExpression::evaluateVariable (const std::string& variableName,
                                               dolfin::Array<double>& values,
                                               const dolfin::Array<double>& x) const
    {
        auto variable = variables_.find (variableName);

        if (variable == variables_.end ())
        {
            dolfin::error ("Cannot find variable \"%s\" in VariableExpression map", variableName.c_str ());
        }

        (variable -> second) -> eval (values, x);
    }
}
