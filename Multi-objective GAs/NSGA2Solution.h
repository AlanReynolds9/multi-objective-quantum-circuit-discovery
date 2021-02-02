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

#ifndef NSGA2_SOLUTION
#define NSGA2_SOLUTION

#include <iostream>

template <class Solution>
class NSGA2Solution
{
  // Solution wrapper for use with NSGA II. This splits those algorithm specific bits - whether the solution is worthy
  // of being in the store, fitness, crowding distances - from those problem specific bits like weights used to
  // calculate objectives, while wrapping those problem specific bits required by the algorithm. Also indicates what is
  // required by the solution class provided by the user.

  // (Should we inherit specific algorithm solutions from a general WrappedSolution class?)
  // (Should the NSGA2Solution class inherit from the problem specific Solution class, hence giving immediate access to
  // the various wrapped methods. (We would lose the fact that this class effectively documents the requirements of a
  // Solution class, but gain brevity!)

public:
  // The templatization in the following function works, but creates a one-to-many relationship, which isn't quite what
  // I was after.
  template <class Sol> friend void crossover(NSGA2Solution<Sol>& lhs, NSGA2Solution<Sol>& rhs);

  // Function to enable sorting and efficient counting of duplicates when calculating population entropy, etc.
  template <class Sol> friend bool sortBefore(const NSGA2Solution<Sol>& lhs, const NSGA2Solution<Sol>& rhs);
  template <class Sol> friend bool lessCrowded(const NSGA2Solution<Sol>& lhs, const NSGA2Solution<Sol>& rhs);

public:
  using Problem = typename Solution::Problem;

  // Wrapped member functions
  NSGA2Solution(const Problem& problem);  // Replaces default constructor - all solutions need to know the problem.
  NSGA2Solution(const NSGA2Solution<Solution>& rhs);
  bool operator==(const NSGA2Solution<Solution>& rhs) const;  // To ensure only unique solutions get into the store.
  bool operator!=(const NSGA2Solution<Solution>& rhs) const;

  void random();
  void mutate(double mutateProb);

  void evaluateObjectives();
  double objective(unsigned obj) const;
  bool dominates(const NSGA2Solution<Solution>& rhs) const;  // Used by the store only.
  bool beats(const NSGA2Solution<Solution>& rhs) const;  // Used by the algorithm. (Beats => dominates.)

  void output(std::ostream& out) const;
  void outputQuality(std::ostream& out) const;

  // NSGA II specific member functions
  int rank() const;
  void setRank(int rank);
  double crowdingDistance() const;
  void setCrowdingDistance(double distance);
  bool worthy() const;

  const Solution& solution() const;  // (Temporary.)

private:
  Solution solution_;         // (Feels like this should perhaps be a pointer??)
  int rank_;
  double crowding_distance;
  mutable bool worthy_;
};

template <class Solution>
std::ostream& operator<<(std::ostream& out, const NSGA2Solution<Solution>& solution);

//----------------------------------------------------------------------------------------------------------------------

// Wrapped functions

template <class Sol>
void crossover(NSGA2Solution<Sol>& lhs, NSGA2Solution<Sol>& rhs)
{
  crossover(lhs.solution_, rhs.solution_);

  // As this has created new solutions, mark them as worthy until it is found that they are dominated.
  // (Necessary as the default copy constructor is used to create children from parents. Might be slightly
  // more efficient to write the copy constructor, but not pleasant.)
  lhs.worthy_ = true;
  rhs.worthy_ = true;
}

template <class Sol>
bool sortBefore(const NSGA2Solution<Sol>& lhs, const NSGA2Solution<Sol>& rhs)
{
  // To enable sorting and efficient counting of duplicates when calculating entropy, etc.
  return sortBefore(lhs.solution_, rhs.solution_);
}

template <class Sol>
bool lessCrowded(const NSGA2Solution<Sol>& lhs, const NSGA2Solution<Sol>& rhs)
{
  return lhs.crowding_distance > rhs.crowding_distance;
}

template <class Sol>
bool fitterThan(const NSGA2Solution<Sol>& lhs, const NSGA2Solution<Sol>& rhs)  // NEW
{
  // Fitness in NSGA2 should be minimized
  if (lhs.rank() != rhs.rank())
  {
    return lhs.rank() < rhs.rank();
  }

  return lessCrowded(lhs, rhs);
}

template <class Solution>
NSGA2Solution<Solution>::NSGA2Solution(const Problem& problem) :
solution_(problem),
rank_(-1),  // An unnatural value - to aid in tracking bugs.
crowding_distance(0),
worthy_(true)
{
}

template <class Solution>
NSGA2Solution<Solution>::NSGA2Solution(const NSGA2Solution<Solution>& rhs) :
solution_(rhs.solution_.clone()),
rank_(-1),  // An unnatural value - to aid in tracking bugs.
crowding_distance(0),
worthy_(true)
{
  // My MOMH code uses this copy constructor solely to create child solutions from parents. It does not necessarily copy
  // 'solution_' from the parent. It may make sense to not copy data generated during evaluation but to only copy the
  // 'genes'. Alternatively, one might give the child a reference to its parent to speed up evaluation. This explains
  // the use of 'clone()' above.
  // (I would prefer write a clone() function for the Solution types associated with each algorithm and then use this in
  // breeder.h and elsewhere. However, this is made trickier by the fact that these Solution types usually don't contain
  // any reference to the Problem being worked on. With no default constructor either, getting the child solution
  // constructed proves problematic.)
}

template <class Solution>
bool NSGA2Solution<Solution>::operator==(const NSGA2Solution<Solution>& rhs) const
{
  return solution_ == rhs.solution_;
}

template <class Solution>
bool NSGA2Solution<Solution>::operator!=(const NSGA2Solution<Solution>& rhs) const
{
  return !operator==(rhs);
}

template <class Solution>
void NSGA2Solution<Solution>::random()
{
  solution_.random();
  worthy_ = true;  // Any new solution is worthy until dominated.
  // (Having to mark solutions as worthy at multiple points in the code seems like a recipe for trouble, much like
  // having to remember to delete items created on the heap.)
}

template <class Solution>
void NSGA2Solution<Solution>::mutate(double mutateProb)
{
  solution_.mutate(mutateProb);
  worthy_ = true;  // Any new solution is worthy until dominated
}

template <class Solution>
void NSGA2Solution<Solution>::evaluateObjectives()
{
  solution_.evaluateObjectives();
}

template <class Solution>
double NSGA2Solution<Solution>::objective(unsigned objNum) const
{
  return solution_.objective(objNum);
}

template <class Solution>
bool NSGA2Solution<Solution>::dominates(const NSGA2Solution<Solution>& rhs) const
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
bool NSGA2Solution<Solution>::beats(const NSGA2Solution<Solution>& rhs) const
{
  // Used by the main GA algorithm. Beats => dominates.
  if (solution_.beats(rhs.solution_))
  {
    rhs.worthy_ = false;
    return true;
  }

  return false;
}

template <class Solution>
void NSGA2Solution<Solution>::output(std::ostream& out) const
{
  solution_.output(out);
}

template <class Solution>
void NSGA2Solution<Solution>::outputQuality(std::ostream& out) const
{
  solution_.outputQuality(out);
}

template <class Solution>
std::ostream& operator<<(std::ostream& out, const NSGA2Solution<Solution>& solution)
{
  solution.output(out);
  return out;
}

// Algorithm specific functions

template <class Solution>
int NSGA2Solution<Solution>::rank() const
{
  return rank_;
}

template <class Solution>
void NSGA2Solution<Solution>::setRank(int rank)
{
  rank_ = rank;
}

template <class Solution>
double NSGA2Solution<Solution>::crowdingDistance() const
{
  return crowding_distance;
}

template <class Solution>
void NSGA2Solution<Solution>::setCrowdingDistance(double distance)
{
  crowding_distance = distance;
}

template <class Solution>
bool NSGA2Solution<Solution>::worthy() const
{
  return worthy_;
}

template <class Solution>  // Temporary
const Solution& NSGA2Solution<Solution>::solution() const
{
  return solution_;
}

#endif
