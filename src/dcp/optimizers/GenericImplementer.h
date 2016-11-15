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

#ifndef SRC_OPTIMIZERS_GENERICIMPLEMENTER_H_INCLUDE_GUARD
#define SRC_OPTIMIZERS_GENERICIMPLEMENTER_H_INCLUDE_GUARD

#include <dcp/problems/GenericEquationSystem.h>
#include <functional>

namespace dcp
{
    /*! \class GenericImplementer GenericImplementer.h
     *  \brief Abstract class that defines the generic API of the classes that implement the specific methods needed by 
     *  the backtracking algorithm. It is template-ized over the type of the control variable
     *  
     */
    template <class T_ControlVariable_>
        class GenericImplementer
        {
            // ---------------------------------------------------------------------------------------------//

            public:
                /************************* TYPEDEFS ************************/
                typedef T_ControlVariable_        T_ControlVariable;
                typedef std::function<void (dcp::GenericEquationSystem&, const T_ControlVariable&)> Updater;
                typedef std::function<T_ControlVariable (const T_ControlVariable&)> SearchDirectionComputer;


                /************************* CONSTRUCTORS ********************/
                //! Constructor [1]
                /*!
                 *  \param updater the updater to be used to update the primal-adjoint system during the backtracking
                 *  algorithm. Must be copy-constructible
                 */
                GenericImplementer (const Updater& updater);

                //! Constructor [2]
                /*!
                 *  \param updater the updater to be used to update the primal-adjoint system during the backtracking
                 *  algorithm. Must be copy-constructible
                 *  \param searchDirectionComputer the object used to compute the search direction. Must be
                 *  copy-constructible
                 */
                GenericImplementer (const Updater& updater, const SearchDirectionComputer& searchDirectionComputer);


                /************************* DESTRUCTOR ***********************/
                //! Destructor
                virtual ~GenericImplementer () {};


                /********************** METHODS ***********************/
                //! Set the dot product to be used
                /*!
                 *  \param dotProductForm the form to be used when computing the dot product 
                 */
                virtual void setDotProduct (const dolfin::Form& dotProductForm);

                //! Set the way the search direction is computed on every loop iteration during minimization
                /*!
                 *  \param searchDirectionComputer the object to be used to compute the search direction
                 */
                virtual void setSearchDirectionComputer (const SearchDirectionComputer& searchDirectionComputer);

                //! Updater for the backtracking algorithm
                /*! 
                 *  This method updates the equations system representing the primal and adjoint problems by using the
                 *  protected member \c updater_
                 *
                 *  \param systems the set of systems to be solved
                 *  \param control the current value of the control function
                 */
                virtual void update (const std::vector<std::shared_ptr<dcp::GenericEquationSystem> > systems,
                                     const T_ControlVariable& control);

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
                                    const std::string& solveType) = 0;

                //! Compute the dot product of the given functions
                /*!
                 *  This actually wraps a call to \c dotProduct_.compute()
                 *  \param left first function
                 *  \param right second function
                 */
                virtual double computeDotProduct (const T_ControlVariable& left,
                                                  const T_ControlVariable& right);

                //! Compute the norm of the given function
                /*!
                 *  This actually wraps a call to \c dotProduct_.norm()
                 *  \param function the function
                 */
                virtual double computeNorm (const T_ControlVariable& function);

                //! Compute the search direction for the backtracking method
                /*!
                 *  This actually wraps a call to \c searchDirectionComputer_
                 *  \param gradient the functional gradient
                 */
                virtual T_ControlVariable computeSearchDirection (const T_ControlVariable& gradient);


                /********************** VARIABLES ***********************/
                //! the problem parameters
                dolfin::Parameters parameters;

                // ---------------------------------------------------------------------------------------------//

            protected:
                //! The updater used to update the system representing primal and ajoint problems
                Updater updater_;

                //! The object used to compute the search direction for every loop iteration.
                /*!
                 *  By default, it is set equal to an object of type dcp::GradientSearchDirection default-constructed,
                 *  but it can be changed using the function \c setSearchDirectionComputer()
                 */
                SearchDirectionComputer searchDirectionComputer_;

                //! The form that will be used to compute the dot product between the gradient and the search direction. 
                /*! 
                 *  The default value is on object of type \c dcp::DotProduct default-constructed, which will try to
                 *  determine the right form to use by checking the geometrical dimensions of the input objects. 
                 *  However, sometimes it may be useful to have a user-defined object to compute the dot product.
                 *  To do so, use the function \c setDotProduct
                 */
                dcp::DotProduct dotProduct_;

                // ---------------------------------------------------------------------------------------------//

            private:

        };



    // ============================================================================================== //
    // ==================================== IMPLEMENTATION ========================================== //
    // ============================================================================================== //

    /******************* CONSTRUCTORS *******************/
    template <class T_ControlVariable>
        GenericImplementer<T_ControlVariable>::GenericImplementer (const Updater& updater) :
            parameters ("generic_implementer_parameters"),
            updater_ (updater),
            searchDirectionComputer_ (dcp::GradientSearchDirection ()),
            dotProduct_ ()
        {
            dolfin::log (dolfin::DBG, "GenericImplementer object created");
        }



    template <class T_ControlVariable>
        GenericImplementer<T_ControlVariable>::GenericImplementer 
                (const Updater& updater, 
                 const SearchDirectionComputer& searchDirectionComputer) :
            parameters ("generic_implementer_parameters"),
            updater_ (updater),
            searchDirectionComputer_ (searchDirectionComputer),
            dotProduct_ ()
        {
            dolfin::log (dolfin::DBG, "GenericImplementer object created");
        }



    /********************** METHODS ***********************/
    template <class T_ControlVariable>
        void GenericImplementer<T_ControlVariable>::setDotProduct (const dolfin::Form& dotProductForm)
        {
            dotProduct_.setDotProductComputer (dotProductForm);
        }



    template <class T_ControlVariable>
        void GenericImplementer<T_ControlVariable>::setSearchDirectionComputer 
            (const SearchDirectionComputer& searchDirectionComputer)
        {
            searchDirectionComputer_ = searchDirectionComputer;
        }



    template <class T_ControlVariable>
        void GenericImplementer<T_ControlVariable>::update 
            (const std::vector<std::shared_ptr<dcp::GenericEquationSystem> > systems,
             const T_ControlVariable& control)
        {
            updater_ (*(systems[0]), control);
        }



    template <class T_ControlVariable>
        double GenericImplementer<T_ControlVariable>::computeDotProduct (const T_ControlVariable& left,
                                                                         const T_ControlVariable& right)
    {
        return dotProduct_.compute (left, right);
    }



    template <class T_ControlVariable>
        double GenericImplementer<T_ControlVariable>::computeNorm (const T_ControlVariable& function)
    {
        return dotProduct_.norm (function);
    }



    template <class T_ControlVariable>
        T_ControlVariable GenericImplementer<T_ControlVariable>::computeSearchDirection 
            (const T_ControlVariable& gradient)
    {
        return searchDirectionComputer_ (gradient);
    }
}
#endif


