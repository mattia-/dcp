/* 
 *  Copyright (C) 2014, Mattia Tamellini, mattia.tamellini@gmail.com
 * 
 *  This file is part of the DCP library
 *   
 *   The DCP library is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   The DCP library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with the DCP library.  If not, see <http://www.gnu.org/licenses/>. 
 */ 

#ifndef SRC_DIFFERENTIAL_PROBLEMS_COMPOSITEDIFFERENTIALPROBLEM_H_INCLUDE_GUARD
#define SRC_DIFFERENTIAL_PROBLEMS_COMPOSITEDIFFERENTIALPROBLEM_H_INCLUDE_GUARD

#include <differential_problems/AbstractProblem.h>
#include <differential_problems/LinearProblem.h>
#include <differential_problems/NonlinearProblem.h>
#include <map>
#include <tuple>
#include <memory>
#include <iostream>
#include <utility>
#include <vector>
#include <functional>

namespace dcp
{
    /*! \class CompositeProblem CompositeProblem.h
     *  \brief Class for multi-variable and multi-equation differential problem
     *  
     *  The class contains a \c std::map that associate a problem with its identifying name
     *  and a vector that stores the problem names in the order they should be solved.
     *  The aforementioned map associates a \c std::string to a pointer to 
     *  \c dcp::AbstractProblem.
     *  The class also offers the possibility to link a problem's coefficient to another problem's
     *  solution through the use of a 
     *  <tt> std::map<std::tuple <std::string, std::string, std::string>, std::pair <std::string, int>> </tt>
     */
    class CompositeProblem
    {
        // ---------------------------------------------------------------------------------------------//  

        public:
            /******************* CONSTRUCTORS ******************/
            //! Default constructor 
            CompositeProblem ();
            
            //! Copy constructor. Deleted, since it makes no sense for unique_ptr's
            CompositeProblem (const CompositeProblem& rhs) = delete;
            
            //! Move constructor
            CompositeProblem (CompositeProblem&& rhs) = default;
            
            
            /******************* DESTRUCTOR *******************/
            //! Default destructor
            ~CompositeProblem () = default;
            

            /******************** OPERATORS *********************/
            //! Copy operator. Deleted since it makes no sense for unique_ptr's
            CompositeProblem& operator= (const CompositeProblem& rhs) = delete;
            
            //! Move operator
            CompositeProblem& operator= (CompositeProblem&& rhs) = default;
            

            /******************** METHODS *********************/
            //! Get composite problem size, i.e. the number of problems stored
            std::size_t size ();
            
            //! Add problem to the map of problems to be solved [1]
            /*!
             *  The parameters are:
             *  \param problemName the problem name
             *  \param problem a const reference to an \c AbstractProblem. 
             *  The class will make a copy of the input problem calling the method \c clone().
             *  The problem's name is inserted at the end of \c solveOrder_
             */
            void addProblem (const std::string& problemName, AbstractProblem& problem);
            
            //! Add problem to the map of problems to be solved [2]
            /*!
             *  The call from main is something like:
             *  \code
             *  std::unique_ptr<dcp::AbstractProblem> foo = 
             *      new dcp::LinearProblem<bilinear_form_type, linear_form_type> (mesh, V));
             *  comp_diff_p.addProblem ("bar", &foo);
             *  \endcode
             *  The parameters are:
             *  \param problemName the problem name
             *  \param problem a unique pointer to a \c dcp::AbstractProblem. 
             *  The class will take ownership of the concrete problem setting the input pointer to \c nullptr.
             *  For example, referring to the snippet of code above, the class will contain (and have full ownership of)
             *  the problem *foo and will set foo to nullptr, so that it cannot be used to modify the class private
             *  members later on.
             *  The problem's name is inserted at the end of \c solveOrder_
             */
            void addProblem (const std::string& problemName, 
                             std::unique_ptr<AbstractProblem>& problem);
            
            //! Remove problem with the given name from the private member variables. A warning is issued if the name is 
            //! not found
            void removeProblem (const std::string& problemName);
            
            //! Set solve order of the problems
            /*!
             *  \param solveOrder a \c std::vector<std::string> in which the problems' names are ordered as one 
             *  wishes the stored problem to be ordered. No check is performed either on problems' names contained in 
             *  the input vector or on vectors' size. This means that, for example, the same problem can be 
             *  insterted more than once, if needed
             */
            void reorderProblems (const std::vector<std::string>& solveOrder);
            
            //! Adds link between problems' coefficient and solution []
            /*!
             *  This function will add to the stored map \c problemsLinks_ a \c std::pair created on the input arguments
             *  and perform the actual linking calling \c linkProblems
             *  \param linkFrom identifies the problem whose parameter (passed as second argument to the function) 
             *  should be linked with the solution of the problem identified by the fourth parameter (\c linkTo)
             *  \param linkedCoefficientName identifies the coefficient to be linked with said solution
             *  \param linkedCoefficientType identifies the type of the coefficient, and will be passed to the function
             *  \c setCoefficient (see \c dcp::AbstractProblem documentation for more details)
             *  \param linkTo identifies the problem whose solution is linked to the parameter in the problem
             *  identified by the second and the first arguments respectively. No check is performed on the
             *  existence of such problem
             *  \param forceRelinking boolean value (default \c false). If the pair (problem name - coefficient) 
             *  identified by the first and the second string in the first argument already appears in the protected 
             *  member variable \c problemsLinks_, it will be relinked using the \c std::pair passed as first argument if 
             *  \c forceRelinking is true, and not relinked if it is false (but issuing a warning in this case)
             */
            void addLink (const std::string& linkFrom, 
                          const std::string& linkedCoefficientName,
                          const std::string& linkedCoefficientType, 
                          const std::string& linkTo,
                          const bool& forceRelinking = false);

            //! Adds link between problems' coefficient and solution [2]
            /*!
             *  This function will add to the stored map \c problemsLinks_ a \c std::pair created on the input arguments
             *  and perform the actual linking calling \c linkProblems
             *  \param linkFrom identifies the problem whose parameter (passed as second argument to the function) 
             *  should be linked with the solution of the problem identified by the fourth parameter (\c linkTo)
             *  \param linkedCoefficientName identifies the coefficient to be linked with said solution
             *  \param linkedCoefficientType identifies the type of the coefficient, and will be passed to the function
             *  \c setCoefficient (see \c dcp::AbstractProblem documentation for more details)
             *  \param linkTo identifies the problem whose solution is linked to the parameter in the problem
             *  identified by the second and the first arguments respectively. No check is performed on the
             *  existence of such problem
             *  \param linkToComponent identifies the component of the solution of the problem indentified by \c linkTo
             *  that should be used as coefficient in the problem identified by \c linkFrom
             *  \param forceRelinking boolean value (default \c false). If the pair (problem name, coefficient) 
             *  identified by the first and the second string in the first argument already appears in the protected 
             *  member variable \c problemsLinks_, it will be relinked using the \c std::pair passed as first argument if 
             *  \c forceRelinking is true, and not relinked if it is false (but issuing a warning in this case)
             */
            void addLink (const std::string& linkFrom, 
                          const std::string& linkedCoefficientName,
                          const std::string& linkedCoefficientType, 
                          const std::string& linkTo,
                          const int& linkToComponent,
                          const bool& forceRelinking = false);
            
            //! Access problem with given name [1] (read only)
            /*!
             *  \param name name of the problem to be accessed. If the name is not found, the function prints an
             *  error message and throws an exception through the function \c dolfin::error()
             *  \return a reference to the problem
             */
            const dcp::AbstractProblem& operator[] (const std::string& name) const;
            
            //! Access problem with given name [2] (read and write)
            /*!
             *  \param name name of the problem to be accessed. If the name is not found, the function prints an
             *  error message and throws an exception through the function \c dolfin::error()
             *  \return a reference to the problem
             */
            dcp::AbstractProblem& operator[] (const std::string& name);
            
            //! Access problem with given position in vector \c solveOrder_ [1] (read only)
            /*!
             *  \param position position of the problem to be accessed in the private member vector \c solveOrder_. 
             *  If \c position is greater than vector size, the function prints an error message and throws an exception 
             *  through the function \c dolfin::error()
             *  \return a reference to the problem
             */
            const dcp::AbstractProblem& operator[] (const std::size_t& position) const;
            
            //! Access problem with given position in vector \c solveOrder_ [2] (read and write)
            /*!
             *  \param position position of the problem to be accessed in the private member vector \c solveOrder_. 
             *  If \c position is greater than vector size, the function prints an error message and throws an exception 
             *  through the function \c dolfin::error()
             *  \return a reference to the problem
             */
            dcp::AbstractProblem& operator[] (const std::size_t& position);
            
            //! Prints information on the problems: names list (in solution order) and links information.
            //! It uses \c dolfin::cout stream
            void print ();
            
            //! Solve all the problems in the order specified by the private member \c solveOrder_
            /*!
             *  \param forceRelinking a boolean flag which, if set to \c true, overrides the current value of protected 
             *  member variable needsLinksScanning_. Default value is \c false
             */
            void solve (const bool& forceRelinking = false);
            
            //! Solve the problem corresponding to the name given [1]
            /*!
             *  \param problemName a string identifying the problem to be solved. If no problem with that name
             *  is found, a warning is issued
             *  \param forceRelinking a boolean flag which, if set to \c true, overrides the current value of protected 
             *  member variable needsLinksScanning_. Default value is \c false
             */
            void solve (const std::string& problemName, const bool& forceRelinking = false);
            
            //! Solve the problem corresponding to the name given [2]
            /*!
             *  This method is provided only to allow calls like
             *  \code
             *  solve ("foo_problem");
             *  \endcode
             *  In this case, the compiler would in fact otherwise call \c solve \c (const \c bool&) which is the 
             *  best-matching implicit conversion for a parameter of type \c const \c char*. Using this method,
             *  the version of \c solve that takes a \c std::string is called as expected.
             */
            void solve (const char* problemName, const bool& forceRelinking = false);
            
            //! Access solution of the problem identified by given name
            /*!
             *  \param problemName name of the problem to be accessed. If the name is not found, the function prints an
             *  error message and throws an exception through the function \c dolfin::error()
             *  \return a reference to the problems' solution
             */
            const dolfin::Function& solution (const std::string& problemName) const;
            
        // ---------------------------------------------------------------------------------------------//  

        protected:
            //! Performs the actual linking between problems
            /*!
             *  Being a protected member, this method is just called from library functions. 
             *  \param link a \c std::pair of the type contained by the protected member \c problemsLinks_
             */
            void linkProblems (const std::pair <
                                                std::tuple <std::string, std::string, std::string>, 
                                                std::pair  <std::string, int>
                                               >& link);
            
            //! The stored problems
            std::map <std::string, std::unique_ptr <dcp::AbstractProblem>> storedProblems_;

            //! The solution order of the problems
            std::vector <std::string> solveOrder_;
            
            //! The map of links between problems. 
            /*!
             *  It is stored as a <tt> map <(string, string, string), (string, int)> </tt> where:
             *  \li the first \c string contains the name of the problem whose coefficient should be linked against 
             *  some other problem's solution
             *  \li the second \c string contains the type of such coefficient, in a form that can be passed to
             *  \c dcp::AbstractProblem::setCoefficients
             *  \li the third \c string contains the name of said coefficient in the problem
             *  \li the fourth \c string contains the name of the problem whose solution should be used to set the 
             *  coefficient identified by the first three strings
             *  \li finally, the last field, which is an \c int, defines which component of the solution of the problem
             *  identified by the fourth \c string should be used in the link. If the whole solution should be used, \c -1
             *  is used as a placeholder
             *  A map guarantees that no tuple (problem, coefficient type, coefficient name) 
             *  is linked twice against possibly different problems
             */
            std::map <std::tuple <std::string, std::string, std::string>, std::pair <std::string, int>> problemsLinks_;
            
        // ---------------------------------------------------------------------------------------------//  

        private:

    };
}

#endif
