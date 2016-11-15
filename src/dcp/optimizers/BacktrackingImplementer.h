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

#ifndef SRC_OPTIMIZERS_BACKTRACKINGIMPLEMENTER_H_INCLUDE_GUARD
#define SRC_OPTIMIZERS_BACKTRACKINGIMPLEMENTER_H_INCLUDE_GUARD

#include <dcp/optimizers/GenericImplementer.h>
#include <dolfin/log/dolfin_log.h>

namespace dcp
{
    /*! \class BacktrackingImplementer BacktrackingImplementer.h
     *  \brief Class that implements the specific methods needed by the backtracking algorithm in the time-independent
     *  case. Derives from \c BacktrackingImplementer .
     *  
     */
    template <class T_ControlVariable>
        class BacktrackingImplementer : public dcp::GenericImplementer<T_ControlVariable>
        {
            // ---------------------------------------------------------------------------------------------//

            public:
                /******************* CONSTRUCTORS *******************/
                //! Constructor [1]
                /*!
                 *  \param updater the updater to be used to update the primal-adjoint system during the backtracking
                 *  algorithm. Must be copy-constructible
                 *
                 *  The constructors also sets the following parameters:
                 *      - \c "primal_problem_name" the name of the object representing the primal problem in the system
                 *        passed to \c apply().
                 *        Default value: "primal"
                 *      - \c "adjoint_problem_name" the name of the object representing the adjoint problem in the system
                 *        passed to \c apply().
                 *        Default value: "adjoint"
                 */
                BacktrackingImplementer 
                    (const typename dcp::GenericImplementer<T_ControlVariable>::Updater& updater);

                //! Constructor [2]
                /*!
                 *  \param updater the updater to be used to update the primal-adjoint system during the backtracking
                 *  algorithm. Must be copy-constructible
                 *  \param searchDirectionComputer the object used to compute the search direction. Must be
                 *  copy-constructible
                 *
                 *  The constructors also sets the following parameters:
                 *      - \c "primal_problem_name" the name of the object representing the primal problem in the system
                 *        passed to \c apply().
                 *        Default value: "primal"
                 *      - \c "adjoint_problem_name" the name of the object representing the adjoint problem in the system
                 *        passed to \c apply().
                 *        Default value: "adjoint"
                 */
                BacktrackingImplementer 
                    (const typename dcp::GenericImplementer<T_ControlVariable>::Updater& updater, 
                     const typename dcp::GenericImplementer<T_ControlVariable>::SearchDirectionComputer& 
                            searchDirectionComputer);


                /************************* DESTRUCTOR ********************/
                //! Destructor
                virtual ~BacktrackingImplementer () {};


                /********************** METHODS ***********************/
                //! Solve the equation systems representing the primal and the adjoint problem. 
                /*!
                 *  \param systems the set of systems to be solved
                 *  \param solveType the type of solve requested; possible values in this class: 
                 *  \li \c all
                 *  \li \c primal
                 *  \li \c adjoint
                 *  with obvious meaning
                 */
                virtual void solve (const std::vector<std::shared_ptr<dcp::GenericEquationSystem> > systems,
                                    const std::string& solveType) override;

                // ---------------------------------------------------------------------------------------------//

            protected:

                // ---------------------------------------------------------------------------------------------//

            private:

        };



    // ============================================================================================== //
    // ==================================== IMPLEMENTATION ========================================== //
    // ============================================================================================== //

    /******************* CONSTRUCTORS *******************/
    template <class T_ControlVariable>
        BacktrackingImplementer<T_ControlVariable>::BacktrackingImplementer
                (const typename dcp::GenericImplementer<T_ControlVariable>::Updater& updater) :
            GenericImplementer<T_ControlVariable> (updater)
        {
            dolfin::begin (dolfin::DBG, "Creating BacktrackingImplementer object...");

            GenericImplementer<T_ControlVariable>::parameters.add ("primal_problem_name", "primal");
            GenericImplementer<T_ControlVariable>::parameters.add ("adjoint_problem_name", "adjoint");

            dolfin::end (); // Creating BacktrackingImplementer object

            dolfin::log (dolfin::DBG, "BacktrackingImplementer object created");
        }



    template <class T_ControlVariable>
        BacktrackingImplementer<T_ControlVariable>::BacktrackingImplementer 
                (const typename dcp::GenericImplementer<T_ControlVariable>::Updater& updater, 
                 const typename dcp::GenericImplementer<T_ControlVariable>::SearchDirectionComputer& 
                        searchDirectionComputer) :
            GenericImplementer<T_ControlVariable> (updater, searchDirectionComputer)
        {
            dolfin::begin (dolfin::DBG, "Creating BacktrackingImplementer object...");

            GenericImplementer<T_ControlVariable>::parameters.add ("primal_problem_name", "primal");
            GenericImplementer<T_ControlVariable>::parameters.add ("adjoint_problem_name", "adjoint");

            dolfin::end (); // Creating BacktrackingImplementer object

            dolfin::log (dolfin::DBG, "BacktrackingImplementer object created");
        }



    /********************** METHODS ***********************/
    template <class T_ControlVariable>
        void BacktrackingImplementer<T_ControlVariable>::solve 
            (const std::vector<std::shared_ptr<dcp::GenericEquationSystem> > systems,
             const std::string& solveType)
        {
            if (solveType == "all")
            {
                systems[0]->solve ();
            }
            else if (solveType == "primal")
            {
                std::string name = GenericImplementer<T_ControlVariable>::parameters["primal_problem_name"];
                (*(systems[0]))[name].solve ();
            }
            else if (solveType == "adjoint")
            {
                std::string name = GenericImplementer<T_ControlVariable>::parameters["adjoint_problem_name"];
                (*(systems[0]))[name].solve ();
            }
            else
            {
                dolfin::dolfin_error ("dcp: BacktrackingImplementer.cpp",
                                      "solve",
                                      "Unknown solve type \"%s\"", 
                                      solveType.c_str ());
            }
        }
}

#endif

