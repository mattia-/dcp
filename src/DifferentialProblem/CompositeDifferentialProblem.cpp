#include <DifferentialProblem/CompositeDifferentialProblem.hpp>
#include <utility>
#include <tuple>
#include <dolfin.h>
#include <algorithm>

namespace control_problem
{

	/******************* METHODS *******************/
	void CompositeDifferentialProblem::addProblem (const std::string& problemName, 
	                                               std::unique_ptr<AbstractDifferentialProblem>& problem)
	{
		// insert problem into storedProblems_ taking ownership
		auto result = storedProblems_.insert (std::make_pair (problemName, std::move (problem)));
		
		// if problem was not inserted in list, issue a warning; else, add it also to the vector problemsOrder_
		if (result.second == false)
		{
			dolfin::warning ("Problem %s already exist in composite differential problem", problemName.c_str ());
		}
		else
		{
			problemsOrder_.emplace_back (problemName);
		}
	}



	void CompositeDifferentialProblem::removeProblem (const std::string& problemName)
	{
		// delete problem from storedProblems_
		auto result = storedProblems_.erase (problemName);
		
		// if problem was not inserted in list, issue a warning; else, add it also to the vector problemsOrder_
		// remember that erase returns the number of elements removed, which in the case of a map is at most 1
		if (result < 1) 
		{
			dolfin::warning ("Problem %s was not removed from composite differential problem", problemName.c_str ());
		}
		else
		{
			// 1)
			// remove problemName from problemsOrder_. Remember that we cannot be sure that such name is 
			// unique in vector and that problems are not in any kind of order.
			// std::find will return an iterator to the element if found and vector::end () if not found
			auto problemPosition = find (problemsOrder_.begin (), problemsOrder_.end (), problemName);
			while (problemPosition != problemsOrder_.end ())
			{
				problemsOrder_.erase (problemPosition);
				problemPosition = find (problemsOrder_.begin (), problemsOrder_.end (), problemName);
			}
			
			// 2)
			// remove problemName from linkedProblems_.
			// We use an important statement from the c++ standard: 
			// When erasing from a map, iterators, pointers and references referring to elements removed by the 
			// function are invalidated. All other iterators, pointers and references keep their validity.
			auto linksIterator = linkedProblems_.begin (); // iterator pointing to the first element of the set
			while (linksIterator != linkedProblems_.end ())
			{
				// delete element if problemName appears either as first or as third element in the tuple
				if (std::get<0> (linksIterator->first) == problemName || linksIterator->second == problemName)
				{
					auto auxIterator = linksIterator; // this will be used for the call to function erase
					++linksIterator;   // iterator incremented to point to next element. 
					                            // This is performed before erasing element, 
					                            // so that increment is still valid
					
					linkedProblems_.erase (auxIterator);
				}
				else
				{
					++linksIterator;
				}
			}
		}
	}



	void CompositeDifferentialProblem::reorderProblems (const std::vector<std::string>& problemsOrder)
	{
		problemsOrder_ = problemsOrder;
	}



	void CompositeDifferentialProblem::linkProblems (const std::string& linkFrom, 
	                                                 const std::string& linkedCoefficientName,
	                                                 const std::string& linkedCoefficientType, 
	                                                 const std::string& linkTo,
	                                                 const bool& forceRelinking = false)
	{
		// create pair containing the link information passed as input arguments.
		// auto keyword used in place of std::pair <std::tuple <std::string, std::string, std::string>, std::string> 
		// to enhance readability
		auto link = std::make_pair (std::make_tuple (linkFrom, linkedCoefficientName, linkedCoefficientType), linkTo);

		// search for map key in linkedProblems_. 
		// remember that the key (i.e. link.first) is an std::tuple<std::string, std::string, std::string>
		// auto keyword used to enhance readability in place of 
		// std::map <std::tuple <std::string, std::string, std::string>, std::string>::iterator 
		auto linkPosition = linkedProblems_.find (link.first);

		if (linkPosition == linkedProblems_.end ()) // if key not found in map, insert link
		{
			linkedProblems_.insert (link);
		}
		else if (forceRelinking == true) // if key found in map but forceRelinking set to true, erase 
			// current link and insert the new one
		{
			dolfin::cout << "In composite control problem: erasing link between problems:" << dolfin::endl;
			dolfin::cout << "\t(" 
				<< std::get<0> (linkPosition->first) 
				<< ", " 
				<< std::get<1> (linkPosition->first) 
				<< ", " 
				<< std::get<2> (linkPosition->first) 
				<< ") -> " 
				<< linkPosition->second 
				<< dolfin::endl;

			linkedProblems_.erase (linkPosition);

			dolfin::cout << "and inserting link: " << dolfin::endl;
			dolfin::cout << "\t(" 
				<< std::get<0> (link.first)
				<< ", " 
				<< std::get<1> (link.first)
				<< ", " 
				<< std::get<2> (link.first)
				<< ") -> " 
				<< link.second 
				<< dolfin::endl;

			linkedProblems_.insert (link);
		}
		else
		{
			dolfin::cout << "Link ("
				<< std::get<0> (link.first)
				<< ", " 
				<< std::get<1> (link.first)
				<< ", " 
				<< std::get<2> (link.first)
				<< ") -> " 
				<< link.second 
				<<") not added because key is already present in stored map" 
				<< dolfin::endl;
		}
	}



	const control_problem::AbstractDifferentialProblem& 
	CompositeDifferentialProblem::problem (const std::string& problemName) const
	{
		auto problemIterator = storedProblems_.find (problemName);
		if (problemIterator == storedProblems_.end ())
		{
			dolfin::error ("Problem %s not found in stored problems map", problemName.c_str ());
		}
		return *(problemIterator->second);
	}



	control_problem::AbstractDifferentialProblem& CompositeDifferentialProblem::problem (const std::string& problemName)
	{
		auto problemIterator = storedProblems_.find (problemName);
		if (problemIterator == storedProblems_.end ())
		{
			dolfin::error ("Problem %s not found in stored problems map", problemName.c_str ());
		}
		return *(problemIterator->second);
	}



	void CompositeDifferentialProblem::print ()
	{
		dolfin::cout << "Problems' solve order:" << dolfin::endl;
		for (auto i : problemsOrder_)
		{
			dolfin::cout << "\t" << i << dolfin::endl; 
		}
		dolfin::cout << dolfin::endl;

		dolfin::cout << "Problems' links:" << dolfin::endl;
		for (auto &i : linkedProblems_)
		{
			dolfin::cout << "("
				<< std::get<0> (i.first)
				<< ", " 
				<< std::get<1> (i.first)
				<< ", " 
				<< std::get<2> (i.first)
				<< ") -> " 
				<< i.second 
				<< dolfin::endl;
		}
	}



	void CompositeDifferentialProblem::solve ()
	{
		// this function iterates over problemsOrder_ and calls solve (problemName) for each problem, thus delegating
		// to the latter function the task of performing the actual parameters setting and solving
		for (auto problem : problemsOrder_)
		{
			solve (problem);
		}
	}



	void CompositeDifferentialProblem::solve (const std::string& problemName)
	{
		dolfin::begin ("Solving problem %s...", problemName.c_str ());
		{	
			// get problem with given name from map. Variable problemIterator will be a
			// std::map <std::string, std::unique_ptr <control_problem::AbstractDifferentialProblem>::iterator
			auto problemIterator = storedProblems_.find (problemName);

			if (problemIterator == storedProblems_.end ())
			{
				dolfin::warning ("Problem %s not found in stored problems map. Aborting...", problemName.c_str ());
				return;
			}

			control_problem::AbstractDifferentialProblem& problem = *(problemIterator->second);

			dolfin::log (dolfin::PROGRESS, "Scanning problems' links...");

			// 1)
			// loop over linkedProblems_. Remember it is a map. Elements in it are order according to the default
			// lexicographical ordering
			auto linksIterator = linkedProblems_.begin ();
			while (linksIterator != linkedProblems_.end () && std::get<0> (linksIterator->first) < problemName)
			{
				if (std::get<0> (linksIterator->first) == problemName)
				{
					// check if target problem of the link exists
					auto targetProblemIterator = storedProblems_.find (linksIterator->second);
					if (targetProblemIterator == storedProblems_.end ())
					{
						dolfin::warning ("Cannot link problem %s to problem %s. No such problem found in \
						                 stored problems map. Aborting...",
						                 problemName.c_str (),
						                 linksIterator->second.c_str ());
						return;
					}

					control_problem::AbstractDifferentialProblem& targetProblem = *(targetProblemIterator->second);

					dolfin::log (dolfin::PROGRESS, 
					             "Linked coefficient %s of type %s to solution of problem %s",
					             (std::get<1> (linksIterator->first)).c_str (),
					             (std::get<2> (linksIterator->first)).c_str (),
					             (linksIterator->second).c_str ());

					problem.setCoefficient (std::get<1> (linksIterator->first), 
					                        boost::shared_ptr<dolfin::Function> (new dolfin::Function (targetProblem.solution ())),
					                        std::get<2> (linksIterator->first));
				}
				++linksIterator;
			}

			// 2)
			// solve problem
			problem.solve ();
		}
		dolfin::end ();
	}
}
