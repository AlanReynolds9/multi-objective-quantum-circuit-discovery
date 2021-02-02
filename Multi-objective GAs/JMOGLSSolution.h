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

#ifndef JMOGLS_SOLUTION
#define JMOGLS_SOLUTION

#include <iostream>

template <class Solution>
class JMOGLSSolution
{
  // Solution wrapper for use with Jasziewicz form of genetic local search. This splits those algorithm specific bits - 
  // whether the solution is worthy of being in the store, fitness - from those problem specific bits like weights used
  // to calculate objectives, while wrapping those problem specific bits required by the algorithm. Also indicates what
  // is required by the solution class provided by the user.

  // (Should we inherit specific algorithm solutions from a general WrappedSolution class?)

public:
  // (The templatization below works, but creates a one-to-many relationship, which isn't quite what I was after.)
  template <class Sol> friend void crossover(JMOGLSSolution<Sol>& lhs, JMOGLSSolution<Sol>& rhs);

  // Function to enable sorting and efficient counting of duplicates when calculating population entropy, etc.
  template <class Sol> friend bool sortBefore(const JMOGLSSolution<Sol>& lhs, const JMOGLSSolution<Sol>& rhs);

  // The oldest solution is removed from the population.
  template <class Sol> friend bool olderThan(const JMOGLSSolution<Sol>& lhs, const JMOGLSSolution<Sol>& rhs);

public:
  using Problem = typename Solution::Problem;

  // Wrapped member functions
  JMOGLSSolution(const Problem& problem);  // Replaces default constructor - all solutions need to know the problem.
  JMOGLSSolution(const JMOGLSSolution<Solution>& rhs);
  bool operator==(const JMOGLSSolution<Solution>& rhs) const;  // To ensure only unique solutions get into the store.
  bool operator!=(const JMOGLSSolution<Solution>& rhs) const;

  void random();
  void mutate(double mutateProb);
  void randomMove();

  void evaluateObjectives();
  double objective(unsigned obj) const;
  bool dominates(const JMOGLSSolution<Solution>& rhs) const;  // Used only by the store.

  void output(std::ostream& out) const;
  void outputQuality(std::ostream& out) const;

  // MOGLS2 specific member functions
  double fitness() const;
  void setFitness(double fitness);

  bool worthy() const;  // (Addition to the store is handled by the algorithm object - this function returns false to...
                        // prevent the survival object from re-adding.) (See definition.)
private:
  static int num_solutions;
  int solution_num;

  Solution solution_;  // (Feels like this should perhaps be a pointer??)
  double fitness_;
};

template <class Solution>
std::ostream& operator<<(std::ostream& out, const JMOGLSSolution<Solution>& solution);

//----------------------------------------------------------------------------------------------------------------------

// Wrapped functions

template <class Sol>
void crossover(JMOGLSSolution<Sol>& lhs, JMOGLSSolution<Sol>& rhs)
{
  crossover(lhs.solution_, rhs.solution_);
}

template <class Sol>
bool sortBefore(const JMOGLSSolution<Sol>& lhs, const JMOGLSSolution<Sol>& rhs)
{
  // To enable sorting and efficient counting of duplicates when calculating entropy, etc.
  return sortBefore(lhs.solution_, rhs.solution_);
}

template <class Solution>
JMOGLSSolution<Solution>::JMOGLSSolution(const Problem& problem) :
solution_num(num_solutions++),
solution_(problem)
{
}

template <class Solution>
JMOGLSSolution<Solution>::JMOGLSSolution(const JMOGLSSolution<Solution>& rhs) :
solution_num(num_solutions++),
solution_(rhs.solution_.clone())
{
  // My MOMH code uses this copy constructor solely to create child solutions from parents. It does not
  // necessarily copy 'solution_' from the parent. It may make sense to not copy data generated during
  // evaluation but to only copy the 'genes'. Alternatively, one might give the child a reference to its
  // parent to speed up evaluation. This explains the use of 'clone()' above.
}

template <class Solution>
bool JMOGLSSolution<Solution>::operator==(const JMOGLSSolution<Solution>& rhs) const
{
  return solution_ == rhs.solution_;
}

template <class Solution>
bool JMOGLSSolution<Solution>::operator!=(const JMOGLSSolution<Solution>& rhs) const
{
  return !operator==(rhs);
}

template <class Solution>
void JMOGLSSolution<Solution>::random()
{
  solution_.random();
}

template <class Solution>
void JMOGLSSolution<Solution>::mutate(double mutateProb)
{
  solution_.mutate(mutateProb);
}

template <class Solution>
void JMOGLSSolution<Solution>::randomMove()
{
  solution_.randomMove();
}

template <class Solution>
void JMOGLSSolution<Solution>::evaluateObjectives()
{
  solution_.evaluateObjectives();
}

template <class Solution>
double JMOGLSSolution<Solution>::objective(unsigned objNum) const
{
  return solution_.objective(objNum);
}

template <class Solution>
bool JMOGLSSolution<Solution>::dominates(const JMOGLSSolution<Solution>& rhs) const
{
  return solution_.dominates(rhs.solution_);
}

template <class Solution>
void JMOGLSSolution<Solution>::output(std::ostream& out) const
{
  solution_.output(out);
}

template <class Solution>
void JMOGLSSolution<Solution>::outputQuality(std::ostream& out) const
{
  solution_.outputQuality(out);
}

template <class Solution>
std::ostream& operator<<(std::ostream& out, const JMOGLSSolution<Solution>& solution)
{
  solution.output(out);
  return out;
}

// Algorithm specific functions

template <class Sol>
bool olderThan(const JMOGLSSolution<Sol>& lhs, const JMOGLSSolution<Sol>& rhs)
{
  return lhs.solution_num < rhs.solution_num;   // Lower solution number implies older
}

template <class Solution>
double JMOGLSSolution<Solution>::fitness() const
{
  return fitness_;
}

template <class Solution>
void JMOGLSSolution<Solution>::setFitness(double fitness)
{
  fitness_ = fitness;
}

template <class Solution>
bool JMOGLSSolution<Solution>::worthy() const
{
  // (This algorithm does not use dominance, so it never determines that a solution is unworthy. However, solutions need
  // to be added to the store on creation, rather than just before removal as is the case with the other algorithms.
  // Hence the addition of solutions to the store is handled by the algorithm object. This function returns 'false' to
  // ensure that the survival object does not try to re-add the solutions upon removal.)
  return false;
}

#endif
