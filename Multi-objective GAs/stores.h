// MOMH-lib (Multi-objective metaheuristic library.)
//
// Copyright (c) 2020  Alan Reynolds (alanreynolds9@gmail.com)
//
// This file is part of MOMH-lib.
//
// MOMH-lib is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
// MOMH-lib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with MOQCD2.  If not, see
// <http://www.gnu.org/licenses/>.

//----------------------------------------------------------------------------------------------------------------------

#ifndef STORES_H
#define STORES_H

// Store classes hold the best solutions found by the algorithm. They may eliminate dominated solutions in each
// generation or at the end of the algorithm, or indeed, upon reaching some maximal size. In either case, before output
// (or other use), the store must ensure that it only holds non-dominated solutions.
//
// Those Store classes that only remove non-dominated solutions before output will be referred to JIT (just-in-time)
// stores. Those that are continuously updated will be referred to as CU or Managed stores.
//
// Worthy solutions are added to the store from the Survival object. Unworthy solutions are marked in the Evaluator
// object. Note that this requires the solution to have public functions bool worthy() and void markUnworthy(). If all
// solutions are not to be added to the store, it is important to ensure that the Evaluator object is used, as not all
// Survival objects require it.

#include <rng.h>


void warnPathologicalCase(bool equal);


// A store that doesn't store anything! Use when using a store would be inefficient or complicated
template <class Solution>
class NoStore
{
public:
  void reset();
  void add(std::shared_ptr<const Solution> addition);  // If solutions are only placed in the store when about to be...
  void sortByObjective(int objNum);                    // ...eliminated from the population, then could we not use...
  void output(std::ostream& out);                      // ...unique_ptrs in all the stores?
  void outputQuality(std::ostream& out);
};

template <class Solution>
std::ostream& operator<<(std::ostream& out, NoStore<Solution>& store);


// A store that can be used regardless of which dominance relation is used, provided that the relation is transitive.
// This means that it can be used with modified dominance relations, but may not be the most efficient when used with
// the standard notion of dominance. JIT means that dominated solutions are only eliminated immediately before printing.
template <class Solution>
class BasicJITStore
{
public:
  void reset();
  void add(std::shared_ptr<const Solution> addition);
  void tidy();  // Now public so that the store can be resorted if the objectives change. (But surely we should only...
                // ...keep those solutions that are best according to the TRUE objectives.)
  void sortByObjective(int objNum);
  void output(std::ostream& out);   // Not const, as store may need to eliminate dominated solutions.
  void outputQuality(std::ostream& out);

private:
  std::vector<std::shared_ptr<const Solution> > member_;
  bool warned_ = false;
};

template <class Solution>
std::ostream& operator<<(std::ostream& out, BasicJITStore<Solution>& store);


// A store that can be used regardless of which dominance relation is used, provided that the relation is transitive,
// that ensures that dominated solutions are never present
template <class Solution>
class BasicManagedStore
{
public:
  void reset();
  void add(std::shared_ptr<const Solution> addition);
  void tidy();  // Added so that the store can be resorted if the objectives change
  void sortByObjective(int objNum);
  void output(std::ostream& out);  // Non-const, as the store is sorted prior to output
  void outputQuality(std::ostream& out);

  // Two functions added to allow MOGLS2 access to some randomly selected elite solutions
  bool empty() const;
  std::shared_ptr<const Solution> randomElement() const;

private:
  std::vector<std::shared_ptr<const Solution> > member_;
  bool warned_ = false;
};

template <class Solution>
std::ostream& operator<<(std::ostream& out, BasicManagedStore<Solution>& store);

//----------------------------------------------------------------------------------------------------------------------

void warnPathologicalCase(bool equal)
{
  // Either two solutions are equal, yet one dominates the other, or two unequal solutions dominate each other. In
  // either case, we warn the user.
  using std::endl;
  using std::cout;

  if (equal)
  {
    cout << "WARNING: Candidate for entering the store dominates a 'guard' solution in the store that" << endl;
    cout << "         appears identical to the candidate. Possibilities include:" << endl << endl;
    cout << "         1. A logically bad dominance relation, where two identical solutions with" << endl;
    cout << "            identical objective values dominate each other." << endl << endl;
    cout << "         2. Operator==() for your solution class does not truly test for equality." << endl;
    cout << "            This may be the case if your solution involves floating point numbers" << endl;
    cout << "            and you consider solutions equal if these numbers match to within some" << endl;
    cout << "            tolerance." << endl << endl;
    cout << "         3. The identical solutions are evaluated in slightly different ways, leading to" << endl;
    cout << "            slightly different objective values." << endl << endl;
    cout << "         The Store class has been fixed to handle this pathological case. (It used to" << endl;
    cout << "         lead to the solutions disappearance from the store.) However, I suggest putting" << endl;
    cout << "         in some serious thought as to whether this might cause issues with other parts" << endl;
    cout << "         of the algorithm." << endl << endl;
  }
  else
  {
    cout << "WARNING: Candidate for entering the store both dominates and is dominated by a" << endl;
    cout << "         (different) guard solution. This suggests that you are using a non-standard" << endl;
    cout << "         dominance relation that does not satisfy the antisymmetry requirement." << endl;
    cout << endl;
  }
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution>
void NoStore<Solution>::reset()
{
}


template <class Solution>
void NoStore<Solution>::add(std::shared_ptr<const Solution> addition)
{
}


template <class Solution>
void NoStore<Solution>::sortByObjective(int objNum)
{
}


template <class Solution>
void NoStore<Solution>::output(std::ostream& out)
{
  out << "No store used." << std::endl;
}


template <class Solution>
void NoStore<Solution>::outputQuality(std::ostream& out)
{
  out << "No store used." << std::endl;
}


template <class Solution>
std::ostream& operator<<(std::ostream& out, NoStore<Solution>& store)
{
  store.output(out);
  return out;
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution>
void BasicJITStore<Solution>::reset()
{
  // If the management of the store accounts for much time, considering reserving and using member_.resize(0) here.
  member_.clear();
}


template <class Solution>
void BasicJITStore<Solution>::add(std::shared_ptr<const Solution> addition)
{
  member_.push_back(addition);
}


template <class Solution>
void BasicJITStore<Solution>::tidy()
{
  // Not the efficient version of Kung et al.

  // Create a temporary vector for storing the non-dominated solutions
  std::vector<std::shared_ptr<const Solution> > nonDominated;
  nonDominated.reserve(member_.size());

  // Compare each solution with those currently in nonDominated.
  for (auto candidate_ptr : member_)
  {
    const auto& candidate = *candidate_ptr;
    bool dominates = false;
    bool dominated = false;

    auto j = nonDominated.begin();
    while (!dominated && j != nonDominated.end())
    {
      const auto& guard = **j;

      // If the candidate has not yet been found to dominate anything, check that it is not dominated.
      // If the candidate is in fact the same solution as an item in the nonDominated list, consider it to be dominated.
      // Remember, a dominated candidate cannot dominate anything in the store of nonDominated solutions, as the
      // dominance relation should be transitive.
      if (!dominates)
      {
        if (guard.dominates(candidate) || guard == candidate)
        {
          dominated = true;
        }
      }

      // If the candidate dominates the guard, remove the guard and note that the candidate is acceptable.
      // It cannot be dominated by anything in the store, as this would mean that the dominating guard must dominate the
      // dominated guard. (Remember, the dominance relation should be transitive.)
      if (candidate.dominates(guard))
      {
        if (dominated)
        {
          // The pathological case where the candidate is equal to or dominated by the guard, but also dominates the
          // guard!!
          if (!warned_)
          {
            warnPathologicalCase(guard == candidate);
            warned_ = true;
          }
        }
        else
        {
          dominates = true;
          j = nonDominated.erase(j);  // (Can we arrage to erase all dominated solution simultaneously? Use worthy_?)
        }
      }
      else
      {
        ++j;
      }
    }

    // If the candidate has survived this long, add it into the vector of nonDominated solutions
    if (!dominated)
    {
      nonDominated.push_back(candidate_ptr);
    }
  }

  // Copy the non-dominated solutions back into member_.  (In modern C++, can we not move them?
  member_ = nonDominated;
  sortByObjective(0);  // (Is this necessary?)
}


template <class Solution>
void BasicJITStore<Solution>::sortByObjective(int objNum)
{
  sort(member_.begin(), member_.end(), CompareByObjective<Solution, std::less>(objNum));
}


template <class Solution>
void BasicJITStore<Solution>::output(std::ostream& out)
{
  tidy();  // (I might like to sort by any objective prior to output. However, this would then resorts by objective 0.)
  for (auto p : member_)
  {
    out << *p << std::endl;
  }
}


template <class Solution>
void BasicJITStore<Solution>::outputQuality(std::ostream& out)
{
  tidy();  // (I might like to sort by any objective prior to output. However, this would then resort by objective 0.)
  for (auto p : member_)
  {
    p->outputQuality(out);
    out << std::endl;
  }
}


template <class Solution>
std::ostream& operator<<(std::ostream& out, BasicJITStore<Solution>& store)
{
  store.output(out);
  return out;
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution>
void BasicManagedStore<Solution>::reset()
{
  // If the management of the store accounts for much time, considering reserving and using member_.resize(0) here.
  member_.clear();
}


template <class Solution>
void BasicManagedStore<Solution>::add(std::shared_ptr<const Solution> addition)
{
  const auto& candidate = *addition;
  bool dominates = false;
  bool dominated = false;

  auto j = member_.begin();
  while (!dominated && j != member_.end())  // If the candidate has been found to be dominated, it cannot dominate...
  {                                         // ...any member of the store.
    const auto& guard = **j;

    // If the candidate has not yet been found to be dominated or dominating, check that it is not dominated by or the
    // same as the guard.
    if (!dominates)
    {
      if (guard.dominates(candidate) || guard == candidate)
      {
        dominated = true;
      }
    }

    // If the candidate dominates the guard, remove the guard and note that the candidate is acceptable.
    if (candidate.dominates(guard))
    {
      if (dominated)
      {
        // The pathological case where the candidate is equal to or dominated by the guard, but also dominates the
        // guard!!
        if (!warned_)
        {
          warnPathologicalCase(guard == candidate);
          warned_ = true;
        }
      }
      else
      {
        dominates = true;
        j = member_.erase(j);  // (Can we arrage to erase all dominated solution simultaneously? Use worthy_?)
      }
    }
    else
    {
      ++j;
    }
  }

  // If the candidate has survived this long, add it into the vector of nonDominated solutions
  if (!dominated)
  {
    member_.push_back(addition);
  }
}


template <class Solution>
void BasicManagedStore<Solution>::tidy()
{
  // Not the efficient version of Kung et al.
  // Used if, for some reason, the dominance relation has changed, eliminating more solutions.

  // Create a temporary vector for storing the non-dominated solutions
  std::vector<std::shared_ptr<const Solution> > nonDominated;
  nonDominated.reserve(member_.size());

  // Compare each solution with those currently in nonDominated.
  for (auto candidate_ptr : member_)
  {
    const auto& candidate = *candidate_ptr;
    bool dominates = false;    // (Repeated code.)
    bool dominated = false;

    auto j = nonDominated.begin();
    while (!dominated && j != nonDominated.end())
    {
      const auto& guard = **j;

      // If the candidate has not yet been found to dominate anything, check that it is not dominated.
      // If the candidate is in fact the same solution as an item in the nonDominated list, consider it to be dominated
      // (Remember, a dominated candidate cannot dominate anything in the store of nonDominated solutions, as the
      // dominance relation should be transitive).
      if (!dominates)
      {
        if (guard.dominates(candidate) || guard == candidate)
        {
          dominated = true;
        }
      }

      // If the candidate dominates the guard, remove the guard and note that the candidate is acceptable.
      // It cannot be dominated by anything in the store, as this would mean that the dominating guard must dominate the
      // dominated guard. (Remember, the dominance relation should be transitive.)
      if (candidate.dominates(guard))
      {
        if (dominated)
        {
          // The pathological case where the candidate is equal to or dominated by the guard, but also dominates the
          // guard!!
          if (!warned_)
          {
            warnPathologicalCase(guard == candidate);
            warned_ = true;
          }
        }
        else
        {
          dominates = true;
          j = nonDominated.erase(j);
        }
      }
      else
      {
        ++j;
      }
    }

    // If the candidate has survived this long, add it into the vector of nonDominated solutions
    if (!dominated)
    {
      nonDominated.push_back(candidate_ptr);
    }
  }

  // Copy the non-dominated solutions back into member_
  member_ = nonDominated;
  sortByObjective(0);  // (Is this necessary?)
}


template <class Solution>
void BasicManagedStore<Solution>::sortByObjective(int objNum)
{
  sort(member_.begin(), member_.end(), CompareByObjective<Solution, std::less>(objNum));
}


template <class Solution>
void BasicManagedStore<Solution>::output(std::ostream& out)
{
  for (auto p : member_)
  {
    out << *p << std::endl;
  }
}


template <class Solution>
void BasicManagedStore<Solution>::outputQuality(std::ostream& out)
{
  for (auto p : member_)
  {
    p->outputQuality(out);
    out << std::endl;
  }
}


template <class Solution>
std::ostream& operator<<(std::ostream& out, BasicManagedStore<Solution>& store)
{
  store.output(out);
  return out;
}


template <class Solution>
bool BasicManagedStore<Solution>::empty() const
{
  return member_.empty();
}


template <class Solution>
std::shared_ptr<const Solution> BasicManagedStore<Solution>::randomElement() const
{
  int r = utils::rand::randInt(0, static_cast<int>(member_.size()));
  return member_[r];
}

//----------------------------------------------------------------------------------------------------------------------

template <class Store, class Iterator>
void addWorthy(Store& store, const Iterator& start, const Iterator& finish)
{
  for (auto i = start; i != finish; ++i)
  {
    if ((*i)->worthy())    // i is an iterator pointing to a std::shared_ptr to a solution
    {
      store.add(*i);
    }
  }
}

#endif
