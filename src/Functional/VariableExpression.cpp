#include <Functional/VariableExpression.hpp>
#include <utility>
#include <memory>

namespace controlproblem
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


    /******************* DESTRUCTOR *******************/
    VariableExpression::~VariableExpression () = default;
    


    /******************* SETTERS *******************/
    void VariableExpression::addVariable (const std::string& variable, 
                                          const boost::shared_ptr <const dolfin::GenericFunction> value,
                                          const bool& forceInsertion)
    {
        dolfin::begin (dolfin::DBG, "Inserting variable in map...");
        if (dolfin::get_log_level () > dolfin::DBG)
        {
            dolfin::end ();
        }
        
        auto result = variables_.insert (std::make_pair (variable, value));
        
        if (result.second == false)
        {
            if (forceInsertion == false)
            {
                dolfin::warning ("Variable \"%s\" already present in map", variable.c_str ());
            }
            else
            {
                dolfin::log (dolfin::DBG, "Updating value of variable \"%s\"...", variable.c_str ());
                variables_.erase (result.first);
                variables_.insert (std::make_pair (variable, value) );
            }
        }
        
        if (dolfin::get_log_level () <= dolfin::DBG)
        {
            dolfin::end ();
        }
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

    

    void VariableExpression::eval (const std::string& variableName,
                                   dolfin::Array<double>& values,
                                   const dolfin::Array<double>& x) const
    {
        dolfin::begin (dolfin::DBG, "Looking for variable \"%s\" in stored map...", variableName.c_str ());
        if (dolfin::get_log_level () > dolfin::DBG)
        {
            dolfin::end ();
        }
        
        auto variable = variables_.find (variableName);
        
        if (variable == variables_.end ())
        {
            dolfin::error ("Cannot find variable \"%s\" in VariableExpression stored map");
        }
        
        (variable -> second) -> eval (values, x);

        if (dolfin::get_log_level () <= dolfin::DBG)
        {
            dolfin::end ();
        }
    }
}
