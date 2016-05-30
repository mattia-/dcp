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
#include <dcp/utils/DotProduct.h>
#include <dolfin/fem/Form.h>
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
             *  \c subiterate_() for more details.
             *  The strings \c storedProblemsSolutionType_ and \c storedProblemsSolveType_ are initialized with
             *  \c "default".
             *  The constructors also sets the following parameters:
             *      - \c "subiterations_tolerance" the tolerance for the convergence check in the subiterations loop.
             *        It will be compared against the sum of the increments' norms divided by the norm of the solution
             *        computed at the last step. Default value: 1e-6
             *      - \c "subiterations_maximum_iterations" the maximum number of iterations allowed in the 
             *        subiterations loop. Default value: 10
             *      - \c "plot_subiteration_solutions": to decide whether to call the \c plotSolution method also during
             *        subiterations. Default: false
             *      - \c "write_subiteration_solutions_to_file": to decide whether to call the \c writeSolutionToFile
             *        method also during subiterations. Default: false
             *      - \c "subiterations_blacklist" a list of problem names whose solution should not be used to check
             *        for the convergence of the subiterations. If empty, all problems in the subiterations range will
             *        be considered. It is overridden by the values in \c "subiterations_whitelist". No check is
             *        performed on the existence of the given names in \c storedProblems_.
             *        The type of the parameters in this nested paremters set does not matter, since
             *        only the keys will be considered. This means that when adding names to this parameter set, any
             *        kind of parameter can be added (even unset parameters), as long as the parameter key matches the 
             *        name of the problem one wishes to blacklist. 
             *        Default value: empty (remember to use parenthesis instead of square brackets to access nested
             *        parameters sets)
             *      - \c "subiterations_whitelist" a list of problem names whose solution should be used to check for
             *        the convergence of the subiterations. If empty, all problems in the subiterations range will be
             *        considered. Overrides values in \c "subiterations_blacklist". No check is performed on the
             *        existence of the given names in \c storedProblems_.
             *        The type of the parameters in this nested paremters set does not matter, since
             *        only the keys will be considered. This means that when adding names to this parameter set, any
             *        kind of parameter can be added (even unset parameters), as long as the parameter key matches the 
             *        name of the problem one wishes to whitelist. 
             *        Default value: empty (remember to use parenthesis instead of square brackets to access nested 
             *        parameters sets)
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
             *  \return a vector containing all the problems names in the order in which they are stored in 
             *  \c solveOrder_ (aka the order in which they are solved
             */
            virtual const std::vector<std::string>& problemsNames () const;
            
            //! Get names of first and last problems on which subiterations will take place
            virtual const std::pair<std::string, std::string>& subiterationsRange () const;
            

            /******************** METHODS *********************/
            //! Add problem to the map of problems to be solved [1]
            /*!
             *  The parameters are:
             *  \param problemName the problem name. Empty names are not allowed.
             *  \param problem a const reference to a \c GenericProblem 
             *  
             *  The class will make a copy of the input problem calling the method \c clone().
             *  The problem's name is inserted at the end of \c solveOrder_
             *
             *  \return \c true if the problem was added, \c false otherwise
             */
            virtual bool addProblem (const std::string& problemName, 
                                     dcp::GenericProblem& problem);
            
            //! Add problem to the map of problems to be solved [2]
            /*!
             *  The parameters are:
             *  \param problemName the problem name. Empty names are not allowed.
             *  \param problem a shared pointer to a \c dcp::GenericProblem. 
             *  
             *  The problem's name is inserted at the end of \c solveOrder_
             *
             *  \return \c true if the problem was added, \c false otherwise
             */
            virtual bool addProblem (const std::string& problemName, 
                                     const std::shared_ptr<dcp::GenericProblem> problem);
            
            //! Add problem to the map of problems to be used to set the initial guesses in the subiteration loop [1]
            /*!
             *  The parameters are:
             *  \param problemName the name of the problem in \c storedProblems_ whose initial guess will be set using
             *  the problem passed as second parameter (no check is performed to determine whether there is a problem
             *  with this name in \c storedProblems_)
             *  \param problem the problem to be used to set the initial guess
             *  
             *  The class will make a copy of the input problem calling the method \c clone().
             *
             *  \return \c true if the problem was added, \c false otherwise
             */
            virtual bool addInitialGuessSetter (const std::string& name, 
                                                dcp::GenericProblem& problem);
            
            //! Add problem to the map of problems to be used to set the initial guesses in the subiteration loop [2]
            /*!
             *  The parameters are:
             *  the problem passed as second parameter (no check is performed to determine whether there is a problem
             *  with this name in \c storedProblems_)
             *  \param problem the problem to be used to set the initial guess
             *  
             *  \return \c true if the problem was added, \c false otherwise
             */
            virtual bool addInitialGuessSetter (const std::string& name, 
                                                const std::shared_ptr<dcp::GenericProblem> problem);
            
            //! Remove problem with the given name from the map where the system problems are stored
            /*! 
             *  Parameters:
             *  \param problemName the name of the problem to be removed
             *
             *  \return \c true if the problem was removed, \c false otherwise
             */
            virtual bool removeProblem (const std::string& problemName);
            
            //! Remove problem with the given name from the map where initial guess setters are stored
            /*! 
             *  Parameters:
             *  \param problemName the name of the problem to be removed
             *
             *  \return \c true if the problem was removed, \c false otherwise
             */
            virtual bool removeInitialGuessSetter (const std::string& name);
            
            //! Set solve order of the problems
            /*!
             *  \param solveOrder a \c std::vector<std::string> in which the problems' names are ordered as one 
             *  wishes the stored problem to be ordered. No check is performed either on problems' names contained in 
             *  the input vector or on vectors' size. The function will, however, check for double entries and empty 
             *  names, since neither of them are allowed. Use \c setSubiterationRange if you want a problems solved more
             *  than once at each timestep
             *
             *  \return \c true if the problems were reordered, \c false otherwise
             */
            virtual bool reorderProblems (const std::vector<std::string>& solveOrder);
            
            //! Adds link between problems' coefficient and solution [1]
            /*!
             *  This function will add to the stored map \c problemsLinks_ a \c std::pair created on the input
             *  arguments.
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
             *
             *  \return \c true if the link was added, \c false otherwise
             */
            virtual bool addLink (const std::string& linkFrom, 
                                  const std::string& linkedCoefficientName,
                                  const std::string& linkedCoefficientType, 
                                  const std::string& linkTo,
                                  const bool& forceRelinking = false);

            //! Adds link between problems' coefficient and solution [2]
            /*!
             *  This function will add to the stored map \c problemsLinks_ a \c std::pair created on the input
             *  arguments.
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
             *
             *  \return \c true if the link was added, \c false otherwise
             */
            virtual bool addLink (const std::string& linkFrom, 
                                  const std::string& linkedCoefficientName,
                                  const std::string& linkedCoefficientType, 
                                  const std::string& linkTo,
                                  const int& linkToComponent,
                                  const bool& forceRelinking = false);
            
            //! Remove link between problems' coefficient and solution
            /*!
             *  Removes the link identified by the input arguments from the protected members.
             *  The input arguments will be used to create an object of \c dcp::GenericProblem::LinkKey to use
             *  to erase the corresponding entry from \c problemsLinks_ . See \c addLink documentation for more details.
             *  The parameter \c linkType has two possible values (\c "system" and \c "initial_guess") which identify
             *  the map from which the link should be removed (respectively, the map of system links used to solve the
             *  system and the map used to set the initial guesses)
             *  
             *  \return \c true if the link was removed, \c false otherwise
             */
            virtual bool removeLink (const std::string& linkFrom, 
                                     const std::string& linkedCoefficientName, 
                                     const std::string& linkedCoefficientType,
                                     const std::string& linkType);
            
            //! Adds link between problems' coefficient and solution for the initial guess setters [1]
            /*!
             *  This function will add to the stored map \c initialGuessesSettersLinks_ a \c std::pair created on the 
             *  input arguments.
             *  See documentation for \c addLink() for more details on the input arguments.
             *
             *  \return \c true if the link was added, \c false otherwise
             */
            virtual bool addInitialGuessSetterLink (const std::string& linkFrom, 
                                                    const std::string& linkedCoefficientName,
                                                    const std::string& linkedCoefficientType, 
                                                    const std::string& linkTo,
                                                    const bool& forceRelinking = false);

            //! Adds link between problems' coefficient and solution [2]
            /*!
             *  This function will add to the stored map \c initialGuessesSettersLinks_ a \c std::pair created on the 
             *  input arguments.
             *  See documentation for \c addLink() for more details on the input arguments.
             *
             *  \return \c true if the link was added, \c false otherwise
             */
            virtual bool addInitialGuessSetterLink (const std::string& linkFrom, 
                                                    const std::string& linkedCoefficientName,
                                                    const std::string& linkedCoefficientType, 
                                                    const std::string& linkTo,
                                                    const int& linkToComponent,
                                                    const bool& forceRelinking = false);
            
            //! Set a new element in \c dotProducts_
            /*!
             *  Adds a new pair in \c dotProducts_ to associate a problem and the \c dcp::DotProduct used to
             *  compute the norm of the increment of its solution.
             *  \param problemName the name of the problem
             *  \param dotProductComputer the dolfin::Form to be used in the \c dcp::DotProduct
             *
             *  \return \c true or \c false indicating the success of the operation
             */
            virtual bool setDotProduct (const std::string& problemName, const dolfin::Form& dotProductComputer);

            //! Access system problem with given name [1] (read only)
            /*!
             *  \param name name of the problem to be accessed. If the name is not found, the function prints an
             *  error message through the function \c dolfin::dolfin_error()
             *  \return a reference to the problem
             */
            virtual const dcp::GenericProblem& operator[] (const std::string& name) const;
            
            //! Access system problem with given name [2] (read and write)
            /*!
             *  \param name name of the problem to be accessed. If the name is not found, the function prints an
             *  error message through the function \c dolfin::dolfin_error()
             *  \return a reference to the problem
             */
            virtual dcp::GenericProblem& operator[] (const std::string& name);
            
            //! Access system problem with given position in vector \c solveOrder_ [1] (read only)
            /*!
             *  \param position position of the problem to be accessed in the private member vector \c solveOrder_. 
             *  If \c position is greater than vector size, the function prints an error message 
             *  through the function \c dolfin::dolfin_error()
             *  \return a reference to the problem
             */
            virtual const dcp::GenericProblem& operator[] (const std::size_t& position) const;
            
            //! Access system problem with given position in vector \c solveOrder_ [2] (read and write)
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
             *  (i.e.: subiterate to the last problem), we use a non-empty string for \c first and an empty string for 
             *  \c last.
             *  No check is performed on the validity of the input strings.
             *  \param first the first problem on which to subiterate (inclusive boundary)
             *  \param last the problem AFTER the last on which to subiterate (exclusive boundary). 
             *  Default value: empty string, which means that we should subiterate up to the last problem.
             */
            virtual void setSubiterationRange (const std::string& first, const std::string& last = "");
            
            //! Prints information on the problems: names list (in solution order) and links information.
            //! It uses \c dolfin::cout stream TODO fix stream and what gets printed
            virtual void print ();
            
            //! Solve all the problems in the order specified by the private member \c solveOrder_
            /*!
             *  This is be performed by calling \c solve() on each problem
             */
            virtual void solve () = 0;
            
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
            
            
            //! Plot the problems' solution
            /*!
             *  \param plotType the plot type desired. It will be passed along to the plot method of each problem, so it
             *  must be a valid string for the problem types stored in the map
             */
            virtual void plotSolution (const std::string& plotType);

            //! Write the problems' solution to file
            /*!
             *  \param writeType the write type desired. It will be passed along to the write method of each problem, so 
             *  it must be a valid string for the problem types stored in the map
             */
            virtual void writeSolutionToFile (const std::string& writeType);
            
            /********************** VARIABLES ***********************/
            //! the system parameters
            dolfin::Parameters parameters;
            
        // ---------------------------------------------------------------------------------------------//  

        protected:
            /******************** METHODS *********************/
            //! Solve the problem corresponding to the name given
            /*!
             *  \param problemName a string identifying the problem to be solved
             */
            virtual void solve_ (const std::string& problemName) = 0;

            //! Add problem to the map of problems to be solved
            /*!
             *  The parameters are:
             *  \param map the map to which the problem should be added
             *  \param problemName the name of the problem 
             *  \param problem a const reference to a \c GenericProblem
             *
             *  \return \c true if the problem was added, \c false otherwise
             */
            virtual bool addProblemToMap_ (std::map<std::string, std::shared_ptr <dcp::GenericProblem>>& map,
                                           const std::string& problemName, 
                                           const std::shared_ptr<dcp::GenericProblem> problem);
            
            //! Remove problem with the given name from the private member variables. A warning is issued if the name is 
            //! not found
            /*!
             *  Parameters:
             *  \param map the map from which the problem should be removed
             *  \param the name of the problem to be removed
             *
             *  \return \c true if the problem was removed, \c false otherwise
             */
            virtual bool removeProblemFromMap_ (std::map<std::string, std::shared_ptr <dcp::GenericProblem>>& map,
                                                const std::string& problemName);
            
            //! Add link to given map
            /*!
             *  Parameters:
             *  \param map the map to which the link should be added
             *  \param link the link to add to the map
             *  \param forceRelinking boolean value. If true, overwrite link
             *
             *  \return \c true if the link was added, \c false otherwise
             */
            virtual bool addLinkToMap_ (std::map<LinkKey, LinkValue>& map, const Link& link, const bool& forceRelinking);

            //! Remove link between problems' coefficient and solution
            /*!
             *  Parameters:
             *  \param map the map from which the link should be removed
             *  \param linkKey the key of link to remove
             *  
             *  \return \c true if the link was removed, \c false otherwise
             */
            virtual bool removeLinkFromMap_ (std::map<LinkKey, LinkValue>& map, const LinkKey& linkKey);
            
            //! Performs the actual linking between problems 
            /*!
             *  Parameters:
             *  \param link a \c std::pair of the type contained by the protected member links maps
             *  \param problemsMap the map in which the problem whose coefficient we want to link is stored
             */
            virtual void linkProblems_ (const Link& link,
                                        std::map<std::string, std::shared_ptr <dcp::GenericProblem>>& problemsMap);
            
            //! Subiterate on given problems
            /*!
             *  The function will subiterate on the problems whose names are stored in \c solveOrder_ 
             *  between \c subiterationsBegin (inclusive) and \c subiterationsEnd (exclusive), 
             *  solving each once in the order provided by \c solveOrder_ until convergence. 
             *  Each problem is solved by calling the method \c solve(problemName)
             *
             *  \param subiterationsBegin iterator to the initial problem on which to subiterate
             *  \param subiterationsEnd iterator to the problem after the last problem on which to subiterate
             */
            //  TODO allow choice of convergence check. For now: max of the norms of the increment and max iter
            virtual void subiterate_ (const std::vector<std::string>::const_iterator subiterationsBegin,
                                      const std::vector<std::string>::const_iterator subiterationsEnd);
                
            //! Auxilary function to set initial guesses for subiterations
            /*! 
             *  \param subiterationsBegin iterator to the initial problem on which to subiterate
             *  \param subiterationsEnd iterator to the problem after the last problem on which to subiterate
             *  \param plotSubiterationSolutions boolean flag
             *  \param writeSubiterationSolutions boolean flag
             */
            virtual void setSubiterationInitialGuesses_ 
                (const std::vector<std::string>::const_iterator subiterationsBegin,
                 const std::vector<std::string>::const_iterator subiterationsEnd,
                 const bool plotSubiterationSolutions,
                 const bool writeSubiterationSolutions);

            //! Auxilary function to solve the problems in the subiterations range
            /*! 
             *  \param subiterationsBegin iterator to the initial problem on which to subiterate
             *  \param subiterationsEnd iterator to the problem after the last problem on which to subiterate
             *  \param plotSubiterationSolutions boolean flag
             *  \param writeSubiterationSolutions boolean flag
             *  \param incrementsNorms vector of the norms of the increments
             *  \param sortedConvergenceCheckProblemNames names of the problems to be used for the convergence check
             *  \oaram oldSolutions vector of the old solutions
             *  \param inLoop boolean flag indicating whether we are in the subiterations loop or before it
             */
            virtual void solveSubiterationProblems_ 
                (const std::vector<std::string>::const_iterator subiterationsBegin,
                 const std::vector<std::string>::const_iterator subiterationsEnd,
                 const int& iteration,
                 const bool plotSubiterationSolutions,
                 const bool writeSubiterationSolutions,
                 std::vector<double>& incrementsNorms,
                 const std::vector<std::string>& sortedConvergenceCheckProblemNames,
                 std::vector<dolfin::Function>& oldSolutions, 
                 const bool inLoop);

            //! Auxilary function to set some subiteration loop variables (the ones that require more instructions)
            /*! 
             *  \param subiterationsBegin iterator to the initial problem on which to subiterate
             *  \param subiterationsEnd iterator to the problem after the last problem on which to subiterate
             *  \param sortedConvergenceCheckProblemNames vector to contain the names of the problems to be used in the
             *  convergence check, sorted as the subiteration sequence
             *  \param oldSolutions vector to contain the solution at the previous subiteration
             */
            virtual void setSubiterationLoopVariables_ (const std::vector<std::string>::const_iterator subiterationsBegin,
                                                        const std::vector<std::string>::const_iterator subiterationsEnd,
                                                        std::vector<std::string>& sortedConvergenceCheckProblemNames,
                                                        std::vector<dolfin::Function>& oldSolutions);

            /******************** VARIABLES *********************/
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
            
            //! A map containing pairing between problems names and the \c dcp::DotProduct s used to compute the norm of 
            //! the increment of the solution of such problems.
            /*!
             *  Basically, the funcion \c subiterate_() looks for the name of the problem that is being solved inside
             *  this map. If found, it uses the \c dcp::DotProduct stored as value in the pair identified by such name
             *  to compute the norm of the increment of the solution of the problem. If not found, it uses the 
             *  default \c dcp::DotProduct class. This map is created empty, and can be populated with the function
             *  \c setDotProductComputer
             */
            std::map <std::string, dcp::DotProduct> dotProducts_;

            //! The problems used to set the initial guesses for the subiterations
            /*! It contains pairs <tt> name - problem </tt>, where \c name identifies a problem in \c storedProblems_
             *  and \c problem contains a \c dcp::GenericProblem that will be solved to set the initial guess for the
             *  problem given by \c name in the subiterations loop. If a problem contained in \c storedProblems_ is not
             *  found in \c initialGuessesSetters_ , then the solution stored in the protected member \c solution_ of
             *  the problem itself (aka: the last computed solution) is used as initial guess.
             */
            std::map<std::string, std::shared_ptr <dcp::GenericProblem>> initialGuessesSetters_;
            
            //! The map of links for the initial guesses setters
            /*! See documentation for \c problemsLinks_ for details. This works exactly in the same way. The only
             *  difference is that while \c problemsLinks_ is used for the problems stored in \c storedProblems_ , 
             *  initialGuessesSettersLinks_ is used for the problems in \c initialGuessesSetters_
             */
            std::map<LinkKey, LinkValue> initialGuessesSettersLinks_;

            //! Solve type to be used in the \c solve() method when solving stored problems
            /*! 
             *  Initial value: \c "default"
             */
            std::string storedProblemsSolveType_;
            
            //! Solution type to be used when calling \c solution() on stored problems
            /*! 
             *  Initial value: \c "default"
             */
            std::string storedProblemsSolutionType_;
        // ---------------------------------------------------------------------------------------------//  

        private:

    };
}

#endif
