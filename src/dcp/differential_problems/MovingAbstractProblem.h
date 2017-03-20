/* 
 *  Copyright (C) 2017, Ivan Fumagalli, ivan.fumagalli@polimi.it
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

#ifndef SRC_DIFFERENTIAL_PROBLEMS_MOVINGABSTRACTPROBLEM_H_INCLUDE_GUARD
#define SRC_DIFFERENTIAL_PROBLEMS_MOVINGABSTRACTPROBLEM_H_INCLUDE_GUARD

// TODO forse in realta' questa classe potrebbe ereditare da qualcosa tipo dcp::EquationSystem

#define EVITAPLOT
#define ADDITIONALVECTOR

#include <dolfin.h>
#include <dcp/differential_problems/SubdomainType.h>
#include <fstream>
#include <dcp/differential_problems/MeshManager.h>

namespace dcp
{

    /*! \class MovingAbstractProblem MovingAbstractProblem.h
     *  \brief Abstract base class for moving-domain differential problems. 
     *
     *  It inherits publicly from \c AbstractProblem and it extends
     *  its functionalities to a moving-domain differential problem.
     *  It is an abstract class, it only provides the
     *  basic interface to all differential problems with moving domains.
     */

class MovingAbstractProblem : virtual public dcp::AbstractProblem
{

  public:

    // Do-nothing non-default constructor
    /*!
     * Ctor required by virtual inheritance from AbstractProblem and
     * possible inheritance from the current class:
     * do-nothing because the present class does not have proper members,
     * non-default because the default ctor is implicitly deleted due
     * to AbstractProblem, and the present ctor needs to call AbstractProblem's
     */
    MovingAbstractProblem (std::shared_ptr<dolfin::FunctionSpace> functionSpace);

    // Return solution function (non const version)
    // TODO sarebbe meglio tenere solo quella const di dcp::AbstractProblem;
    //      per il momento questa serve per poter passare *solution_.vector() a dolfin::LinearSolver
    dolfin::Function& solution ();


	  //! Performs operations before moving the mesh
	  /*!
	   *	Purely virtual method, to allow the assembling of part of the system on the previous domain
	   */
    virtual void preassemble () = 0;

    //! Set mesh-and-dofs manager
    virtual void setMeshManager (const MeshManager<dolfin::ALE> & meshManager);

    //! Get mesh-and-dofs manager
    virtual dcp::MeshManager<>  meshManager () const
      { return * meshManager_; }

    //! Update the problem after mesh changing
    /*! After the mesh has changed, the forms, function spaces and functions defined on it have to be updated.
     *  This update is performed via dolfin::adapt
     */
    virtual void adapt () = 0;
    //TODO remove! if not necessary

    //! Clone method
    /*!
     *  \return a pointer to a \c dcp::MovingAbstractProblem containing a copy of the object on 
     *  which it is called. 
     *  Pure virtual: comes as pure virtual also from dcp::AbstractProblem
     */

  protected:

    const MeshManager<dolfin::ALE> * meshManager_;

}; //end of class definition

} //end of namespace
#endif
