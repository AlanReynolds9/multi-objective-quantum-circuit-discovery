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

#ifndef MOGA_SOLUTION
#define MOGA_SOLUTION

#include <iostream>

template <class Solution>
class MOGASolution
{
  // Solution wrapper for use with MOGA. This splits those algorithm specific bits - whether the solution is worthy
  // of being in the store, fitness, crowding distances - from those problem specific bits like weights used to
  // calculate objectives, while wrapping those problem specific bits required by the algorithm. Also indicates what is
  // required by the solution class provided by the user.

  // (Should we inherit specific algorithm solutions from a general WrappedSolution class?)

public:
  // The templatization of the following function works, but creates a one-to-many relationship, which isn't quite what
  // I was after.
  template <class Sol> friend void crossover(MOGASolution<Sol>& lhs, MOGASolution<Sol>& rhs);

  // Function to enable sorting and efficient counting of duplicates when calculating population entropy, etc.
  template <class Sol> friend bool sortBefore(const MOGASolution<Sol>& lhs, const MOGASolution<Sol>& rhs);

public:
  using Problem = typename Solution::Problem;

  // Wrapped member functions
  MOGASolution(const Problem& problem);  // Replaces default constructor - all solutions need to know about the problem.
  MOGASolution(const MOGASolution<Solution>& rhs);
  bool operator==(const MOGASolution<Solution>& rhs) const;   // To ensure only unique solutions get into the store.
  bool operator!=(const MOGASolution<Solution>& rhs) const;

  void random();
  void mutate(double mutateProb);

  void evaluateObjectives();
  double objective(unsigned obj) const;
  bool dominates(const MOGASolution<Solution>& rhs) const;  // Used by the store only.
  bool beats(const MOGASolution<Solution>& rhs) const;  // Used by the algorithm. (Beats => dominates.)

  void output(std::ostream& out) const;
  void outputQuality(std::ostream& out) const;

  // MOGA II specific member functions
  double fitness() const;
  void setFitness(double fitness);
  bool worthy() const;

private:
  Solution solution_;         // (Feels like this should perhaps be a pointer??)
  double fitness_;
  mutable bool worthy_;
};

template <class Solution>
std::ostream& operator<<(std::ostream& out, const MOGASolution<Solution>& solution);

//----------------------------------------------------------------------------------------------------------------------

// Wrapped functions

template <class Sol>
void crossover(MOGASolution<Sol>& lhs, MOGASolution<Sol>& rhs)
{
  crossover(lhs.solution_, rhs.solution_);

  // As this has created new solutions, mark them as worthy until it is found that they are dominated.
  // (Necessary as the default copy constructor is used to create children from parents. Might be slightly
  // more efficient to write the copy constructor, but not pleasant.)
  lhs.worthy_ = true;
  rhs.worthy_ = true;
}

template <class Sol>
bool sortBefore(const MOGASolution<Sol>& lhs, const MOGASolution<Sol>& rhs)
{
  // To enable sorting and efficient counting of duplicates when calculating entropy, etc.
  return sortBefore(lhs.solution_, rhs.solution_);
}

template <class Solution>
MOGASolution<Solution>::MOGASolution(const Problem& problem) :
solution_(problem),
worthy_(true)
{
}

template <class Solution>
MOGASolution<Solution>::MOGASolution(const MOGASolution<Solution>& rhs) :
solution_(rhs.solution_.clone()),
worthy_(true)
{
  // My MOMH code uses this copy constructor solely to create child solutions from parents. It does not
  // necessarily copy 'solution_' from the parent. It may make sense to not copy data generated during
  // evaluation but to only copy the 'genes'. Alternatively, one might give the child a reference to its
  // parent to speed up evaluation. This explains the use of 'clone()' above.
}

template <class Solution>
bool MOGASolution<Solution>::operator==(const MOGASolution<Solution>& rhs) const
{
  return solution_ == rhs.solution_;
}

template <class Solution>
bool MOGASolution<Solution>::operator!=(const MOGASolution<Solution>& rhs) const
{
  return !operator==(rhs);
}

template <class Solution>
void MOGASolution<Solution>::random()
{
  solution_.random();
  worthy_ = true;  // Any new solution is worthy until dominated
}

template <class Solution>
void MOGASolution<Solution>::mutate(double mutateProb)
{
  solution_.mutate(mutateProb);
  worthy_ = true;  // Any new solution is worthy until dominated
}

template <class Solution>
void MOGASolution<Solution>::evaluateObjectives()
{
  solution_.evaluateObjectives();
}

template <class Solution>
double MOGASolution<Solution>::objective(unsigned objNum) const
{
  return solution_.objective(objNum);
}

template <class Solution>
bool MOGASolution<Solution>::dominates(const MOGASolution<Solution>& rhs) const
{
  // Used by the Store only.
  if (solution_.dominates(rhs.solution_))
  {
    rhs.worthy_ = false;
    return true;
  }

  return false;
}

template <class Solution>
bool MOGASolution<Solution>::beats(const MOGASolution<Solution>& rhs) const
{
  // Used by the main GA algorithm.
  if (solution_.beats(rhs.solution_))
  {
    rhs.worthy_ = false;
    return true;
  }

  return false;
}

template <class Solution>
void MOGASolution<Solution>::output(std::ostream& out) const
{
  solution_.output(out);
}

template <class Solution>
void MOGASolution<Solution>::outputQuality(std::ostream& out) const
{
  solution_.outputQuality(out);
}

template <class Solution>
std::ostream& operator<<(std::ostream& out, const MOGASolution<Solution>& solution)
{
  solution.output(out);
  return out;
}

// Algorithm specific functions

template <class Solution>
double MOGASolution<Solution>::fitness() const
{
  return fitness_;
}

template <class Solution>
void MOGASolution<Solution>::setFitness(double fitness)
{
  fitness_ = fitness;
}

template <class Solution>
bool MOGASolution<Solution>::worthy() const
{
  return worthy_;
}

#endif
