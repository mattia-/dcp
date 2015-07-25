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

#ifndef SRC_PROBLEMS_GENERICEQUATIONSYSTEM_H_INCLUDE_GUARD
#define SRC_PROBLEMS_GENERICEQUATIONSYSTEM_H_INCLUDE_GUARD

#include <dcp/problems/GenericProblem.h>
#include <map>
#include <tuple>
#include <memory>
#include <iostream>
#include <utility>
#include <vector>
#include <functional>

namespace dcp
{
    /*! \class GenericEquationSystem GenericEquationSystem.h
     *  \brief Generic class for multi-variable and multi-equation coupled system
     *  
     *  The class contains a \c std::map that associate a problem with its identifying name
     *  and a vector that stores the problem names in the order they should be solved.
     *  The aforementioned map associates a \c std::string to a pointer to 
     *  \c dcp::GenericProblem.
     *  The class also offers the possibility to link a problem's coefficient to another problem's
     *  solution through the use of a 
     *  <tt> std::map<std::tuple <std::string, std::string, std::string>, std::pair <std::string, int>> </tt>.
     */
    class GenericEquationSystem
    {
        // ---------------------------------------------------------------------------------------------//  

        public:
            /******************* TYPEDEFS **********************/
            typedef std::tuple <std::string, std::string, std::string> LinkKey;
            typedef std::pair <std::string, int> LinkValue;
            typedef std::pair <LinkKey, LinkValue> Link;
            

            /******************* CONSTRUCTORS ******************/
            //! Defaul constructor 
            /*!
             *  The strings in \c subiterationsRange_ are initialized with empty strings. 
             *  This way, no subiteration is performed (unless \c setSubiterationRange() is explicitly called), 
             *  since empty name are not allowed for problems. See documentation for \c setSubiterationRange() and
             *  \c subiterate() for more details.
             *  The strings \c solutionType_ and \c solveType_ are initialized with \c "default".
             *  The constructors also sets the following parameters:
             *      - \c "subiterations_tolerance" the tolerance for the convergence check in the subiterations loop.
             *        It will be compared against the sum of the increments' norms. Default value: 1e-6
             *      - \c "subiterations_maximum_iterations" the maximum number of iterations allowed in the 
             *        subiterations loop. Default value: 10
             */
            GenericEquationSystem ();
            
            
            /******************* DESTRUCTOR *******************/
            //! Default destructor
            virtual ~GenericEquationSystem () = default;
            

            /******************** GETTERS *********************/
            //! Get equation system size, i.e. the number of problems stored
            virtual const std::size_t size () const;
            
            //! Get names of problems stored in the protected map \c storedProblems_
            /*!
             *  \return a vector containing all the problems names. Move semantic is 
             *  used in returning the vector, so no unnecessary copy is performed.
             */
            virtual const std::vector<std::string> problemsNames () const;
            
            //! Get names of first and last problems on which subiterations will take place
            virtual const std::pair<std::string, std::string>& subiterationsRange () const;
            

            /******************** METHODS *********************/
            //! Add problem to the map of problems to be solved [1]
            /*!
             *  The parameters are:
             *  \param problemName the problem name. Empty names are not allowed.
             *  \param problem a const reference to an \c GenericProblem 
             *  
             *  The class will make a copy of the input problem calling the method \c clone().
             *  The problem's name is inserted at the end of \c solveOrder_
             */
            virtual void addProblem (const std::string& problemName, dcp::GenericProblem& problem);
            
            //! Add problem to the map of problems to be solved [2]
            /*!
             *  The parameters are:
             *  \param problemName the problem name. Empty names are not allowed.
             *  \param problem a shared pointer to a \c dcp::GenericProblem. 
             *  
             *  The problem's name is inserted at the end of \c solveOrder_
             */
            virtual void addProblem (const std::string& problemName, 
                                     const std::shared_ptr<dcp::GenericProblem> problem);
            
            //! Remove problem with the given name from the private member variables. A warning is issued if the name is 
            //! not found
            virtual void removeProblem (const std::string& problemName);
            
            //! Set solve order of the problems
            /*!
             *  \param solveOrder a \c std::vector<std::string> in which the problems' names are ordered as one 
             *  wishes the stored problem to be ordered. No check is performed either on problems' names contained in 
             *  the input vector or on vectors' size. The function will, however, check for double entries and empty 
             *  names, since neither of them are not allowed
             */
            virtual void reorderProblems (const std::vector<std::string>& solveOrder);
            
            //! Adds link between problems' coefficient and solution [1]
            /*!
             *  This function will add to the stored map \c problemsLinks_ a \c std::pair created on the input arguments
             *  and perform the actual linking calling \c linkProblems
             *  \param linkFrom identifies the problem whose parameter (passed as second argument to the function) 
             *  should be linked with the solution of the problem identified by the fourth parameter (\c linkTo)
             *  \param linkedCoefficientName identifies the coefficient to be linked with said solution
             *  \param linkedCoefficientType identifies the type of the coefficient, and will be passed to the function
             *  \c setCoefficient (see \c dcp::GenericProblem documentation for more details)
             *  \param linkTo identifies the problem whose solution is linked to the parameter in the problem
             *  identified by the second and the first arguments respectively. No check is performed on the
             *  existence of such problem
             *  \param forceRelinking boolean value (default \c false). If the tuple identifying the coefficient already
             *  appears in the protected member variable \c problemsLinks_, it will be relinked using the \c std::pair
             *  passed as first argument if \c forceRelinking is true, and not relinked if it is false (but issuing a
             *  warning in this case)
             */
            virtual void addLink (const std::string& linkFrom, 
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
             *  \c setCoefficient (see \c dcp::GenericProblem documentation for more details)
             *  \param linkTo identifies the problem whose solution is linked to the parameter in the problem
             *  identified by the second and the first arguments respectively. No check is performed on the
             *  existence of such problem
             *  \param linkToComponent identifies the component of the solution of the problem indentified by \c linkTo
             *  that should be used as coefficient in the problem identified by \c linkFrom
             *  \param forceRelinking boolean value (default \c false). If the tuple identifying the coefficient already
             *  appears in the protected member variable \c problemsLinks_, it will be relinked using the \c std::pair
             *  passed as first argument if \c forceRelinking is true, and not relinked if it is false (but issuing a
             *  warning in this case)
             */
            virtual void addLink (const std::string& linkFrom, 
                                  const std::string& linkedCoefficientName,
                                  const std::string& linkedCoefficientType, 
                                  const std::string& linkTo,
                                  const int& linkToComponent,
                                  const bool& forceRelinking = false);
            
            //! Remove link between problems' coefficient and solution
            /*!
             *  Removes the link identified by the input arguments from the protected member \c problemsLinks_ .
             *  The input arguments will be used to create an object of \c dcp::GenericProblem::LinkKey to use
             *  to erase the corresponding entry from \c problemsLinks_ .
             *  
             *  \return \c true if the link was removed, \c false otherwise
             */
            virtual bool removeLink (const std::string& linkFrom, 
                                     const std::string& linkedCoefficientName, 
                                     const std::string& linkedCoefficientType);
            
            //! Access problem with given name [1] (read only)
            /*!
             *  \param name name of the problem to be accessed. If the name is not found, the function prints an
             *  error message through the function \c dolfin::dolfin_error()
             *  \return a reference to the problem
             */
            virtual const dcp::GenericProblem& operator[] (const std::string& name) const;
            
            //! Access problem with given name [2] (read and write)
            /*!
             *  \param name name of the problem to be accessed. If the name is not found, the function prints an
             *  error message through the function \c dolfin::dolfin_error()
             *  \return a reference to the problem
             */
            virtual dcp::GenericProblem& operator[] (const std::string& name);
            
            //! Access problem with given position in vector \c solveOrder_ [1] (read only)
            /*!
             *  \param position position of the problem to be accessed in the private member vector \c solveOrder_. 
             *  If \c position is greater than vector size, the function prints an error message 
             *  through the function \c dolfin::dolfin_error()
             *  \return a reference to the problem
             */
            virtual const dcp::GenericProblem& operator[] (const std::size_t& position) const;
            
            //! Access problem with given position in vector \c solveOrder_ [2] (read and write)
            /*!
             *  \param position position of the problem to be accessed in the private member vector \c solveOrder_. 
             *  If \c position is greater than vector size, the function prints an error message 
             *  through the function \c dolfin::dolfin_error()
             *  \return a reference to the problem
             */
            virtual dcp::GenericProblem& operator[] (const std::size_t& position);
            
            //! Set subiteration range
            /*! Basically, \c subiterate will be called on the range of problems in \c solveOrder_ going from
             *  \c subiterationsRange_.first (inclusive) to \c subiterationsRange_.second (exclusive). 
             *  The two strings in \c subiterationsRange_ are initialized as empty strings, so that if unless
             *  \c setSubiterationRange() is called explicitly no subiteration is performed, since
             *  empty names cannot be used for problems in \c storedProblems_. To indicate the end of the vector
             *  (i.e.: subiterate to the last problem), use a non-empty string for \c first and an empty string for 
             *  \c last.
             *  No check is performed on the validity of the input strings.
             */
            virtual void setSubiterationRange (const std::string& first, const std::string& last);
            
            //! Prints information on the problems: names list (in solution order) and links information.
            //! It uses \c dolfin::cout stream TODO fix stream and what gets printed
            virtual void print ();
            
            //! Solve all the problems in the order specified by the private member \c solveOrder_
            /*!
             *  This is be performed by calling \c solve() on each problem
             */
            virtual void solve () = 0;
            
            //! Solve the problem corresponding to the name given
            /*!
             *  \param problemName a string identifying the problem to be solved. If no problem with that name
             *  is found, a warning is issued
             */
            virtual void solve (const std::string& problemName) = 0;
            
            //! Access solution of the problem identified by given name
            /*!
             *  \param problemName name of the problem to be accessed. If the name is not found, the function prints an
             *  error message through the function \c dolfin::dolfin_error()
             *  \param solutionType the solution type requested. Passed along to \c solution() when called on the 
             *  problem identified by \c problemName. Default value: \c "default" 
             *  
             *  \return a reference to the problems' solution
             */
            virtual const dolfin::Function& solution (const std::string& problemName, 
                                                      const std::string& solutionType = "default") const;
            
            
            /********************** VARIABLES ***********************/
            //! the system parameters
            dolfin::Parameters parameters;
            
        // ---------------------------------------------------------------------------------------------//  

        protected:
            /******************** METHODS *********************/
            //! Subiterate on given problems
            /*!
             *  The function will subiterate on the problems whose names are stored in \c solveOrder_ 
             *  between \c subiterationsBegin (inclusive) and \c subiterationsEnd (exclusive), 
             *  solving each once in the order provided by \c solveOrder_ until convergence. 
             *  Each problem is solved by calling the method \c solve(problemName)
             *  TODO allow choice of convergence check. For now: sum of the norms of the increment and max iter
             */
            virtual void subiterate (std::vector<std::string>::const_iterator subiterationsBegin,
                                     std::vector<std::string>::const_iterator subiterationsEnd);
                
            /******************** VARIABLES *********************/
            //! Performs the actual linking between problems
            /*!
             *  Being a protected member, this method is just called from library functions. 
             *  \param link a \c std::pair of the type contained by the protected member \c problemsLinks_
             */
            virtual void linkProblems (const Link& link);
            
            //! The stored problems
            std::map<std::string, std::shared_ptr <dcp::GenericProblem>> storedProblems_;

            //! The solution order of the problems
            std::vector<std::string> solveOrder_;
            
            //! The map of links between problems. 
            /*!
             *  It is stored as a <tt> map <(string, string, string), (string, int)> </tt> where:
             *  \li the first \c string contains the name of the problem whose coefficient should be linked against 
             *  some other problem's solution
             *  \li the second \c string contains the type of such coefficient, in a form that can be passed to
             *  \c dcp::GenericProblem::setCoefficients
             *  \li the third \c string contains the name of said coefficient in the problem
             *  \li the fourth \c string contains the name of the problem whose solution should be used to set the 
             *  coefficient identified by the first three strings
             *  \li finally, the last field, which is an \c int, defines which component of the solution of the problem
             *  identified by the fourth \c string should be used in the link. If the whole solution should be used, 
             *  \c -1 is used as a placeholder
             *  A map guarantees that no tuple (problem, coefficient type, coefficient name) 
             *  is linked twice against possibly different problems
             */
            std::map<LinkKey, LinkValue> problemsLinks_;
            
            //! The strings specifying the name of the first and last problem on which we must subiterate
            std::pair<std::string, std::string> subiterationsRange_;
            
            //! Solve type to be used in the \c solve() method when solving stored problems
            /*! 
             *  Initial value: \c "default"
             */
            std::string solveType_;
            
            //! Solution type to be used when calling \c solution() on stored problems
            /*! 
             *  Initial value: \c "default"
             */
            std::string solutionType_;
        // ---------------------------------------------------------------------------------------------//  

        private:

    };
}

#endif
