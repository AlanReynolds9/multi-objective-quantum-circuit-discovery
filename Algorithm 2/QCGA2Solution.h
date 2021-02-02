// MOQCD2 (Multi-objective quantum circuit discovery - algorithm 2.)
//
// Copyright (c) 2020  Alan Reynolds (alanreynolds9@gmail.com)
//
// This file is part of MOQCD2.
//
// MOQCD2 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
// MOQCD2 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with MOQCD2.  If not, see
// <http://www.gnu.org/licenses/>.

//----------------------------------------------------------------------------------------------------------------------

#ifndef QCGA2_SOLUTION
#define QCGA2_SOLUTION

#include <iostream>

template <class Solution>
class QCGA2Solution
{
  // Solution wrapper for use with QCGA2. This splits those algorithm specific bits - whether the solution is worthy of
  // being in the store, fitness, crowding distances - from those problem specific bits like weights used to calculate
  // objectives, while wrapping those problem specific bits required by the algorithm. Also indicates what is required
  // by the solution class provided by the user.

  // Essentially the same as NSGA2Solution, but with minor changes.

public:
  // The templating on the function below works, but creates a one-to-many relationship, which isn't quite what I was
  // after.
  template <class Sol> friend void crossover(QCGA2Solution<Sol>& lhs, QCGA2Solution<Sol>& rhs);

  // Enable sorting and efficient counting of duplicates. Used when calculation population entropy, etc.
  template <class Sol> friend bool sortBefore(const QCGA2Solution<Sol>& lhs, const QCGA2Solution<Sol>& rhs);

  template <class Sol> friend bool lessCrowded(const QCGA2Solution<Sol>& lhs, const QCGA2Solution<Sol>& rhs);

public:
  using Problem = typename Solution::Problem;

  // Wrapped member functions
  QCGA2Solution(const Problem& problem);  // Replaces default constructor. All solutions need to know about the problem.
  QCGA2Solution(const QCGA2Solution<Solution>& rhs);
  bool operator==(const QCGA2Solution<Solution>& rhs) const;   // To ensure only unique solutions get into the store.
  bool operator!=(const QCGA2Solution<Solution>& rhs) const;

  void random();
  void mutate(double mutateProb);

  void evaluateObjectives();
  bool bad() const;
  double objective(unsigned obj) const;
  bool dominates(const QCGA2Solution<Solution>& rhs) const;  // Used by the store only.
  bool beats(const QCGA2Solution<Solution>& rhs) const;  // Used by the algorithm. (Beats => dominates.)

  void output(std::ostream& out) const;
  void outputQuality(std::ostream& out) const;  // (We could supply numObj() instead, allowing stores and populations...
                                                // ...to fetch objective values when outputing quality.)
  // QCGA II specific member functions
  int rank() const;
  void setRank(int rank);
  double crowdingDistance() const;
  void setCrowdingDistance(double distance);
  bool worthy() const;

  const Solution& solution() const;  // Temporary?

private:
  Solution solution_;         // (It feels like this should perhaps be a pointer??)
  int rank_;
  double crowding_distance;
  mutable bool worthy_;
};

template <class Solution>
std::ostream& operator<<(std::ostream& out, const QCGA2Solution<Solution>& solution);

//----------------------------------------------------------------------------------------------------------------------

// Wrapped functions

template <class Sol>
void crossover(QCGA2Solution<Sol>& lhs, QCGA2Solution<Sol>& rhs)
{
  // No need to set 'worthy_' - the copy constructor does this when the child is cloned from the parent.
  crossover(lhs.solution_, rhs.solution_);
}

template <class Sol>
bool sortBefore(const QCGA2Solution<Sol>& lhs, const QCGA2Solution<Sol>& rhs)
{
  // To enable sorting and efficient counting of duplicates when calculating entropy, etc.
  return sortBefore(lhs.solution_, rhs.solution_);
}

template <class Sol>
bool lessCrowded(const QCGA2Solution<Sol>& lhs, const QCGA2Solution<Sol>& rhs)
{
  return lhs.crowding_distance > rhs.crowding_distance;
}

template <class Sol>
bool fitterThan(const QCGA2Solution<Sol>& lhs, const QCGA2Solution<Sol>& rhs)  // NEW
{
  // Fitness in QCGA2 should be minimized
  if (lhs.rank() != rhs.rank())
  {
    return lhs.rank() < rhs.rank();
  }

  return lessCrowded(lhs, rhs);
}

template <class Solution>
QCGA2Solution<Solution>::QCGA2Solution(const Problem& problem) :
solution_(problem),
rank_(-1),  // An unnatural value - to aid in tracking bugs.
crowding_distance(0),
worthy_(true)
{
}

template <class Solution>
QCGA2Solution<Solution>::QCGA2Solution(const QCGA2Solution<Solution>& rhs) :
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
bool QCGA2Solution<Solution>::operator==(const QCGA2Solution<Solution>& rhs) const
{
  return solution_ == rhs.solution_;
}

template <class Solution>
bool QCGA2Solution<Solution>::operator!=(const QCGA2Solution<Solution>& rhs) const
{
  return !operator==(rhs);
}

template <class Solution>
void QCGA2Solution<Solution>::random()
{
  // No need to set 'worthy_' - the Solutions constructor does this already.
  solution_.random();
}

template <class Solution>
void QCGA2Solution<Solution>::mutate(double mutateProb)
{
  // No need to set 'worthy_' - the copy constructor does this when the child is cloned from the parent.
  // Note that mutateProb is always 1.1 - i.e. always mutate.
  solution_.mutate();
}

template <class Solution>
void QCGA2Solution<Solution>::evaluateObjectives()
{
  solution_.evaluateObjectives();
}


template <class Solution>
bool QCGA2Solution<Solution>::bad() const
{
  return solution_.bad();
}


template <class Solution>
double QCGA2Solution<Solution>::objective(unsigned objNum) const
{
  return solution_.objective(objNum);
}

template <class Solution>
bool QCGA2Solution<Solution>::dominates(const QCGA2Solution<Solution>& rhs) const
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
bool QCGA2Solution<Solution>::beats(const QCGA2Solution<Solution>& rhs) const
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
void QCGA2Solution<Solution>::output(std::ostream& out) const
{
  solution_.output(out);
}

template <class Solution>
void QCGA2Solution<Solution>::outputQuality(std::ostream& out) const
{
  solution_.outputQuality(out);
}

template <class Solution>
std::ostream& operator<<(std::ostream& out, const QCGA2Solution<Solution>& solution)
{
  solution.output(out);
  return out;
}

// Algorithm specific functions

template <class Solution>
int QCGA2Solution<Solution>::rank() const
{
  return rank_;
}

template <class Solution>
void QCGA2Solution<Solution>::setRank(int rank)
{
  rank_ = rank;
}

template <class Solution>
double QCGA2Solution<Solution>::crowdingDistance() const
{
  return crowding_distance;
}

template <class Solution>
void QCGA2Solution<Solution>::setCrowdingDistance(double distance)
{
  crowding_distance = distance;
}

template <class Solution>
bool QCGA2Solution<Solution>::worthy() const
{
  return worthy_;
}

template <class Solution>  // Temporary?
const Solution& QCGA2Solution<Solution>::solution() const
{
  return solution_;
}

#endif
