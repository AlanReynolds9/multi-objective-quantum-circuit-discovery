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

#ifndef POPULATION_H
#define POPULATION_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <functional>
#include <vector>
#include <algorithm>
#include <memory>
//#include <execution>  // Not yet supported on either Clang or gcc.
#include <future>
//#include <thread>
#include "utils.h"

// Much of the time, when writing Survival classes, etc., using a simple vector is much more straightforward. We should
// consider adding converions from vector to Population and vice-versa. Also a  move constructor might be useful.

template <class Solution>
class Population
{
public:
  // Constructors, copying and merging
  explicit Population(int cap);
  Population(const Population<Solution>& rhs);
  Population(const Population<Solution>& lhs, const Population<Solution>& rhs); // Merge constructor

  Population<Solution>& operator=(const Population<Solution>& rhs);

  // Element access
  using SolutionType = Solution;
  using iterator = typename std::vector<std::shared_ptr<Solution> >::iterator;
  using const_iterator = typename std::vector<std::shared_ptr<Solution> >::const_iterator;
  iterator begin();
  const_iterator cbegin() const;
  iterator end();
  const_iterator cend() const;
  
  std::shared_ptr<Solution>& operator[](int i);  // (Is it possible to use a unique_ptr here?)
  const std::shared_ptr<const Solution> operator[](int i) const;  // (Either use const and &, or neither?)

  // Creation of an initial population.
  void initialize(const typename Solution::Problem& problem, int numSolutions);  // For partial filling of the population
  void initialize(const typename Solution::Problem& problem);  // Filling the whole population

  // Evaluation of the objectives for the entire population.
  void evaluateObjectives();

  // Addition of a solution to a non-full population
  void add(std::shared_ptr<Solution> addition);

  // Emptying a population
  void clear();

  // Population size
  int capacity() const; // (Added for the SPEA2 survival class. Probably a better choice than some of the calls to...
  int size() const;     // ...size() in the older algorithms.)
  bool full() const;

  // Output
  void output(std::ostream& out) const;
  void outputQuality(std::ostream& out) const;
  double entropy() const;

private:
  int capacity_;
  std::vector<std::shared_ptr<Solution> > member_;  // (Is it possible to use a unique_ptr here?)
};  

//----------------------------------------------------------------------------------------------------------------------

// Must include the implementation in the header file in order to be able to use templated classes from another file.

template <class Solution>
Population<Solution>::Population(int cap):
capacity_(cap)
{
  // Creates a population with a fixed capacity
  member_.reserve(cap);
}


template <class Solution>
Population<Solution>::Population(const Population<Solution>& rhs):
capacity_(rhs.capacity_),
member_(rhs.member_)
{
  // This copies pointers to solutions, not the solutions themselves.
}


template <class Solution>
Population<Solution>::Population(const Population<Solution>& lhs, const Population<Solution>& rhs):
capacity_(lhs.capacity_ + rhs.capacity_),
member_(lhs.member_)
{
  // A merge constructor. Constructs a population containing all the elements of lhs and rhs.
  member_.insert(member_.end(), rhs.member_.begin(), rhs.member_.end()); // (Hmm, seems to work if I remove member_...
}                                                                        // ...from the rhs. bit)


template <class Solution>
Population<Solution>& Population<Solution>::operator=(const Population<Solution>& rhs)
{
  // Check for assignment to self
  if (this == &rhs)
  {
    return *this;
  }
  
  // Copy the population of pointers to solutions.
  capacity_ = rhs.capacity_;
  member_ = rhs.member_;
  return *this;
}


template <class Solution>
typename Population<Solution>::iterator Population<Solution>::begin()
{
  return member_.begin();
}


template <class Solution>
typename Population<Solution>::const_iterator Population<Solution>::cbegin() const
{
  return member_.begin();
}


template <class Solution>
typename Population<Solution>::iterator Population<Solution>::end()
{
  return member_.end();
}


template <class Solution>
typename Population<Solution>::const_iterator Population<Solution>::cend() const
{
  return member_.end();
}


// (I would like to be able to remove the first of the access methods. This allows you to update the solutions. We would
// like to be able to store the total fitness of the population in the Population class, but this will become invalid if
// one can change solutions willy nilly.)
template <class Solution>
std::shared_ptr<Solution>& Population<Solution>::operator[](int i)
{
  return member_[i];
}


template <class Solution>
const std::shared_ptr<const Solution> Population<Solution>::operator[](int i) const
{
  return member_[i];
}


template <class Solution>
void Population<Solution>::initialize(const typename Solution::Problem& problem, int numSolutions)
{
  // Fill the population, perhaps only partially, with pointers to randomly generated solutions.
  assert(numSolutions <= capacity_);

  clear();
  for (auto i = 0; i < numSolutions; ++i)
  {
    auto solution(std::make_shared<Solution>(problem));
    solution->random();
    member_.push_back(solution);
   }

  // Evaluate all the solutions in the population.
  evaluateObjectives();
}


template <class Solution>
void Population<Solution>::initialize(const typename Solution::Problem& problem)
{
  // Fill the population with pointers to randomly generated solutions. 
  initialize(problem, capacity_);
}


// There follows a number of different versions of evaluateObjectives(). The first is the serial version which is handy
// for debugging and testing purposes. The second is a too modern version that my compilers can't cope with yet. The
// third is an alternative manual parallelization, that works if each solution takes approximately the same amount of
// time to evaluate. Alas, this is not the case for much of my work. The fourth, uncommented version, uses a pool of
// work that threads grab as they become idle. This is the one I currently use for my experiments.

//template <class Solution>  // Serial version
//void Population<Solution>::evaluateObjectives()
//{
//  for_each(member_.begin(), member_.end(), [](auto& solution){solution->evaluateObjectives();});
//}


//template <class Solution>  // C++17 version - parallel algorithms not yet supported on Clang or gcc.
//void Population<Solution>::evaluateObjectives()
//{
//  using std::execution::par;  // Not yet supported on Clang or gcc.
//  for_each(par, member_.begin(), member_.end(), [](auto& solution){solution->evaluateObjectives();});
//}


//template <class Solution>
//void Population<Solution>::evaluateObjectives()
//{
//  // This version determines the number, n, of hardware threads we can have. It then divides the population into n
//  // approximately equal chunks, hopefully sharing the work out evenly. Each task is set off using std::async.
//  // (Note that using std::async for each individual (unevaluated) solution is slower, due to the overheads associated
//  // with setting up each thread. (This last note is based on testing with old quantum circuit optimization code.)
//  int n = std::thread::hardware_concurrency();
//  std::vector<std::future<void> > futures;
//  futures.reserve(n);
//  auto lambda = [](iterator begin, iterator end)
//                  {for_each(begin, end, [](auto& solution){solution->evaluateObjectives();});};
//  auto policy = std::launch::async;  // Tried using 'deferred' for just the first - didn't make any difference.
//  for (auto i = 0; i < n; ++i)
//  {
//    futures.push_back(std::async(policy, lambda, member_.begin() + i * size() / n,
//                                 member_.begin() + (i + 1) * size() / n));
//  }
//  for (auto& f : futures)
//  {
//    f.get();
//  }
//}


template <class Solution>  // Function object for the manual thread approach to parallelizing evaluateObjectives().
class WorkFinder
{
  using Element = std::shared_ptr<Solution>;
  using Iterator = Element*;
public:
  WorkFinder(std::vector<std::shared_ptr<Solution> >& pool) :
  next_(pool.data()),
  end_(pool.data() + pool.size())
  {
  }

  void operator()()
  {
    while (true)
    {
      auto current = next_++;
      if (current >= end_)
      {
        return;
      }
      (*current)->evaluateObjectives();
    }
  }

private:
  std::atomic<Iterator> next_;
  Iterator end_;
};


// Manual thread approach to parallelizing evaluateObjectives. (A tiny bit faster than using std:async - perhaps.)
template <class Solution>
void Population<Solution>::evaluateObjectives()
{
  WorkFinder<Solution> doMyShare(member_);

  std::vector<std::thread> threads;
  for (auto i = 0; i < std::thread::hardware_concurrency() - 1; ++i)
  {
    std::thread t{std::ref(doMyShare)};
    threads.push_back(std::move(t));
  }
  doMyShare();

  for (auto i = 0; i < std::thread::hardware_concurrency() - 1; ++i)
  {
    threads[i].join();
  }
}

template <class Solution>
void Population<Solution>::add(std::shared_ptr<Solution> addition)
{
  assert(!full());
  member_.push_back(addition);
}

template <class Solution>
void Population<Solution>::clear()
{
  member_.resize(0);    // Avoiding VC++.NET's habit of reducing a vector's capacity to zero on use of 'clear()'.
}

template <class Solution>
int Population<Solution>::capacity() const
{
  return capacity_;
}

template <class Solution>
int Population<Solution>::size() const
{
  return static_cast<int>(member_.size());
}

template <class Solution>
bool Population<Solution>::full() const
{
  return size() == capacity_;
}

template <class Solution>
void Population<Solution>::output(std::ostream& out) const
{
  for (const_iterator i = member_.cbegin(); i != member_.cend(); ++i)  // Code compiles (on the Mac) using just...
  {                                                                    // ...begin() and end().
    out << **i << std::endl;
  }
}

template <class Solution>
void Population<Solution>::outputQuality(std::ostream& out) const
{
  for (const_iterator i = member_.cbegin(); i != member_.cend(); ++i)
  {
    (*i)->outputQuality(out);
    out << std::endl;
  }
}

template <class Solution>
std::ostream& operator<<(std::ostream& out, const Population<Solution>& population)
{
  population.output(out);
  return out;
}

template <class Solution>
double Population<Solution>::entropy() const
{
  // Calculate the entropy of the population.
  // This assumes that suitable equality and sortBefore functions have been provided for the solutions

  // Copy the solutions into a vector that we can sort, and sort it.
  std::vector<std::shared_ptr<Solution> > popCopy(member_);  // auto?
  sort(popCopy.begin(), popCopy.end(), sortBeforeP<std::shared_ptr<Solution> >);

  // Now count the number of times each solution occurs in the population
  // (It is interesting to wonder whether using upper_bound woule be faster than using find_if. Binary search is likely
  // to be more efficient when there is low entropy and less efficient with high entropy.)
  std::vector<int> counts;
  const_iterator i = popCopy.cbegin();
  const_iterator e = popCopy.cend();
  while (i != e)
  {
    // Count the number of solutions that are the same as that pointed to by i. (Actually pointed to by the pointer...
    // ...pointed to by i!)
    const_iterator next = find_if(i, e, PointerNotEqual<std::shared_ptr<Solution> >(*i));  // auto?
    int count = static_cast<int>(next - i);

    // Take a note
    counts.push_back(count);

    i = next;
  }

  double ent = 0;
  for (auto i = 0; i < counts.size(); ++i)
  {
    double prob = static_cast<double>(counts[i]) / size();
    ent -= prob * log(prob) / log(2.0);
  }

  return ent;
}

#endif
