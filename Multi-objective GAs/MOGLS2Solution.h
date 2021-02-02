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

#ifndef MOGLS2_SOLUTION
#define MOGLS2_SOLUTION

#include <iostream>

template <class Solution>
class MOGLS2Solution
{
  // Solution wrapper for use with Ishibuchi and Murata's ammended genetic local search. This splits those algorithm
  // specific bits - whether the solution is worthy of being in the store, fitness - from those problem specific bits
  // like weights used to calculate objectives, while wrapping those problem specific bits required by the algorithm.
  // Also indicates what is required by the solution class provided by the user.

  // (Should we inherit specific algorithm solutions from a general WrappedSolution class?)

public:
  // The templatization in the following function works, but creates a one-to-many relationship, which isn't quite what
  // I was after.
  template <class Sol> friend void crossover(MOGLS2Solution<Sol>& lhs, MOGLS2Solution<Sol>& rhs);

  // Function to enable sorting and efficient counting of duplicates when calculating population entropy, etc.
  template <class Sol> friend bool sortBefore(const MOGLS2Solution<Sol>& lhs, const MOGLS2Solution<Sol>& rhs);

public:
  using Problem = typename Solution::Problem;

  // Wrapped member functions
  MOGLS2Solution(const Problem& problem); // Replaces default constructor - all solutions need to know about the problem
  MOGLS2Solution(const MOGLS2Solution<Solution>& rhs);
  bool operator==(const MOGLS2Solution<Solution>& rhs) const;   // To ensure only unique solutions get into the store.
  bool operator!=(const MOGLS2Solution<Solution>& rhs) const;

  void random();
  void mutate(double mutateProb);
  void randomMove();

  void evaluateObjectives();
  double objective(unsigned obj) const;
  bool dominates(const MOGLS2Solution<Solution>& rhs) const;  // Used only by the store.

  void output(std::ostream& out) const;
  void outputQuality(std::ostream& out) const;

  // MOGLS2 specific member functions
  double fitness() const;
  void setFitness(double fitness);
  double angle() const;
  void setAngle(double angle);

  bool worthy() const;  // (Addition to the store is handled by the algorithm object - this function merely returns
                        // ...false to prevent the survival object from re-adding.) (See function definition.)
private:
  Solution solution_;         // (Feels like this should perhaps be a pointer??)
  double fitness_;
  double angle_;  // Used to determine search direction in local search
};

template <class Solution>
std::ostream& operator<<(std::ostream& out, const MOGLS2Solution<Solution>& solution);

//----------------------------------------------------------------------------------------------------------------------

// Wrapped functions

template <class Sol>
void crossover(MOGLS2Solution<Sol>& lhs, MOGLS2Solution<Sol>& rhs)
{
  crossover(lhs.solution_, rhs.solution_);
}

template <class Sol>
bool sortBefore(const MOGLS2Solution<Sol>& lhs, const MOGLS2Solution<Sol>& rhs)
{
  // To enable sorting and efficient counting of duplicates when calculating entropy, etc.
  return sortBefore(lhs.solution_, rhs.solution_);
}

template <class Solution>
MOGLS2Solution<Solution>::MOGLS2Solution(const Problem& problem) :
solution_(problem),
angle_(0)
{
}

template <class Solution>
MOGLS2Solution<Solution>::MOGLS2Solution(const MOGLS2Solution<Solution>& rhs) :
solution_(rhs.solution_.clone()),
angle_(0)
{
}

template <class Solution>
bool MOGLS2Solution<Solution>::operator==(const MOGLS2Solution<Solution>& rhs) const
{
  return solution_ == rhs.solution_;
}

template <class Solution>
bool MOGLS2Solution<Solution>::operator!=(const MOGLS2Solution<Solution>& rhs) const
{
  return !operator==(rhs);
}

template <class Solution>
void MOGLS2Solution<Solution>::random()
{
  solution_.random();
}

template <class Solution>
void MOGLS2Solution<Solution>::mutate(double mutateProb)
{
  solution_.mutate(mutateProb);
}

template <class Solution>
void MOGLS2Solution<Solution>::randomMove()
{
  solution_.randomMove();
}

template <class Solution>
void MOGLS2Solution<Solution>::evaluateObjectives()
{
  solution_.evaluateObjectives();
}

template <class Solution>
double MOGLS2Solution<Solution>::objective(unsigned objNum) const
{
  return solution_.objective(objNum);
}

template <class Solution>
bool MOGLS2Solution<Solution>::dominates(const MOGLS2Solution<Solution>& rhs) const
{
  return solution_.dominates(rhs.solution_);
}

template <class Solution>
void MOGLS2Solution<Solution>::output(std::ostream& out) const
{
  solution_.output(out);
}

template <class Solution>
void MOGLS2Solution<Solution>::outputQuality(std::ostream& out) const
{
  solution_.outputQuality(out);
}

template <class Solution>
std::ostream& operator<<(std::ostream& out, const MOGLS2Solution<Solution>& solution)
{
  solution.output(out);
  return out;
}

// Algorithm specific functions

template <class Solution>
double MOGLS2Solution<Solution>::fitness() const
{
  return fitness_;
}

template <class Solution>
void MOGLS2Solution<Solution>::setFitness(double fitness)
{
  fitness_ = fitness;
}

template <class Solution>
double MOGLS2Solution<Solution>::angle() const
{
  return angle_;
}

template <class Solution>
void MOGLS2Solution<Solution>::setAngle(double angle)
{
  angle_ = angle;
}

template <class Solution>
bool MOGLS2Solution<Solution>::worthy() const
{
  // The algorithm does not use dominance, so it never determines that a solution is unworthy. However, solutions need
  // to be added to the store on creation, rather than just before removal as is the case with the other algorithms.
  // Hence the addition of solutions to the store is handled by the algorithm object. This function returns 'false' to
  // ensure that the survival object does not try to re-add the solutions upon removal.
  return false;
}

#endif
