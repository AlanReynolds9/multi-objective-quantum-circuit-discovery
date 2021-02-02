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

#ifndef BREEDERS_H
#define BREEDERS_H

#include <array>  // NEW
#include <vector>
#include <algorithm>
#include <memory>
#include <rng.h>

// An alternative approach to the breeder classes might be to use inheritance, creating a `Breeder' abstract base class,
// an `AdjustableBreeder' abstract class that has the methods to adjust mutation and crossover rates, etc. (Could still
// templatize the algorithm class on the Breeder class.)

// (Might there be justification for changing child() so that it returns a Solution by value? Well, perhaps it should be
// cheap to move Solutions, but since we probably don't have any guarantees.)

template <class Solution, class Selector>
class StandardBreeder
{
public:
  explicit StandardBreeder(Selector& selector);
  std::unique_ptr<Solution> child();  // Breeder owns child until it is passed to algorithm.
  void reset();  // Do we need this, or could we just ensure that a new breeder is created in each iteration?

  void setMutationProb(double newProb);
  void setCrossoverProb(double newProb);

private:
  void fill();
  void breed();

private:
  Selector& selector_;

  double mutation_prob;
  double crossover_prob;

  std::array<std::unique_ptr<Solution>, 2> child_;  // Breeder owns the child until it is passed out to the algorithm.
  int child_index{0};
};

template <class Solution, class Selector>
class SingleChildBreeder
{
  // Same as a standard breeder, except that only one child is created from each pair of parents. Use the same approach
  // of filling the bed, breeding and child selection as for the standard breeder.
  // (Could we implement this in terms of a StandardBreeder??)
public:
  explicit SingleChildBreeder(Selector& selector);
  std::unique_ptr<Solution> child();
  void reset();   // Do we need this, or could we just ensure that a new breeder is created in each iteration

  void setMutationProb(double newProb);
  void setCrossoverProb(double newProb);

private:
  void fill();
  void breed();

private:
  Selector& selector_;

  double mutation_prob;
  double crossover_prob;

  // Bed has room for two solutions, but only the first child is ever selected
  std::array<std::unique_ptr<Solution>, 2> child_;
};

template <class Solution, class Selector>
class OrgyBreeder
{
  // Performs all selection at once, filling a mating pool that is typically the same size as the population.
  // (Adding mating restrictions to this would be tricky. Would probably want to create a new class from scratch.)
public:
  explicit OrgyBreeder(Selector& selector);
  std::unique_ptr<Solution> child();
  void reset();   // Do we need this, or could we just ensure that a new breeder is created in each iteration

  void setSize(int newSize);
  void setMutationProb(double newProb);
  void setCrossoverProb(double newProb);

private:
  void fill();
  void breed();

private:
  Selector& selector_;

  int size_;
  double mutation_prob;
  double crossover_prob;

  // We cannot use an std::array here, as the size of the structure is unknown until the user calls setSize. While we
  // could initialize the following in-class, it's probably neater in this case to do it in the constructor.
  std::vector<std::unique_ptr<Solution> > child_;  // Cannot use a std::array here. (Size unknown)
  int child_index;
};

template <class Solution, class Selector>
class GPBreeder
{
  // Variant of the StandardBreeder to be used for genetic programming, where either crossover is performed, or
  // mutation.
public:
  explicit GPBreeder(Selector& selector);
  std::unique_ptr<Solution> child();
  void reset();    // Do we need this, or could we just ensure that a new breeder is created in each iteration

  void setCrossoverProb(double newProb);   // If not crossover, then mutation

private:
  void fill();
  void breed();

private:
  Selector& selector_;

  bool performing_crossover{false};
  double crossover_prob;

  std::array<std::unique_ptr<Solution>, 2> child_;
  int child_index{0};
};

//----------------------------------------------------------------------------------------------------------------------

template <class Solution, class Selector>
StandardBreeder<Solution, Selector>::StandardBreeder(Selector& selector):
selector_(selector)  // (Mutation and crossover probabilities are uninitialized!)
{
}

template <class Solution, class Selector>
std::unique_ptr<Solution> StandardBreeder<Solution, Selector>::child()
{
  // If we are getting the first child, first create the children required.
  if (child_index == 0)
  {
    fill();
    breed();
  }

  // Return the referred to child, updating the child_index to point to the next
  auto oldIndex = child_index;
  child_index = 1 - child_index;

  return std::move(child_[oldIndex]);  // Can't copy a unique_ptr - hence std::move.
}

template <class Solution, class Selector>
void StandardBreeder<Solution, Selector>::reset()
{
  // Resets the breeder, between iterations, to ensure that new children are created the next time 'child' is called.
  child_index = 0;
}

template <class Solution, class Selector>
void StandardBreeder<Solution, Selector>::setMutationProb(double newProb)
{
  // MUST use this function before using the Breeder.
  // (Is it possible to include the setting of mutatation and crossover probabilities to user defined values in the
  // constructor, or is it essential that all Breeder constructors have the same signature?)
  mutation_prob = newProb;
}

template <class Solution, class Selector>
void StandardBreeder<Solution, Selector>::setCrossoverProb(double newProb)
{
  // MUST use this function before using the Breeder. See above.
  crossover_prob = newProb;
}

template <class Solution, class Selector>
void StandardBreeder<Solution, Selector>::fill()
{
  // Creates copies of the selected parent solutions. This, followed by crossover, is somewhat inefficient.
  // Also inefficient to create and then copy.
  child_[0] = std::make_unique<Solution>(*selector_.select());
  child_[1] = std::make_unique<Solution>(*selector_.select());
}

template <class Solution, class Selector>
void StandardBreeder<Solution, Selector>::breed()
{
  if (utils::rand::rand01() < crossover_prob)
  {
    crossover(*child_[0], *child_[1]);
  }
  child_[0]->mutate(mutation_prob);
  child_[1]->mutate(mutation_prob);
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution, class Selector>
SingleChildBreeder<Solution, Selector>::SingleChildBreeder(Selector& selector):
selector_(selector)  // (Mutation and crossover probabilities are uninitialized!)
{
}

template <class Solution, class Selector>
std::unique_ptr<Solution> SingleChildBreeder<Solution, Selector>::child()
{
  // New parents are required whenever we create a child.
  fill();
  breed();

  return std::move(child_[0]);
}

template <class Solution, class Selector>
void SingleChildBreeder<Solution, Selector>::reset()
{
  // Since only a single child is created at a time, there is no need to reset.
}

template <class Solution, class Selector>
void SingleChildBreeder<Solution, Selector>::setMutationProb(double newProb)
{
  mutation_prob = newProb;
}

template <class Solution, class Selector>
void SingleChildBreeder<Solution, Selector>::setCrossoverProb(double newProb)
{
  crossover_prob = newProb;
}

template <class Solution, class Selector>
void SingleChildBreeder<Solution, Selector>::fill()
{
  // Creates copies of the selected parent solutions. This, followed by crossover, is somewhat inefficient.
  // Also inefficient to create and then copy.
  child_[0] = std::make_unique<Solution>(*selector_.select());  // Creates a new solution, using copy constructor.
  child_[1] = std::make_unique<Solution>(*selector_.select());
}

template <class Solution, class Selector>
void SingleChildBreeder<Solution, Selector>::breed()
{
  if (utils::rand::rand01() < crossover_prob)
  {
    crossover(*child_[0], *child_[1]);
  }
  child_[0]->mutate(mutation_prob);
  // No need to mutate the second child, as it is culled.
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution, class Selector>
OrgyBreeder<Solution, Selector>::OrgyBreeder(Selector& selector):
selector_(selector),
size_(0),
child_(0),
child_index(0)  // (Mutation and crossover probabilities are uninitialized!)
{
}

template <class Solution, class Selector>
std::unique_ptr<Solution> OrgyBreeder<Solution, Selector>::child()
{
  // Should only use to get size_ children
  if (child_index == size_)
  {
    std::cerr << "Attempting to acquire too many children from an OrgyBreeder. Aborting." << std::endl;
    exit(EXIT_FAILURE);
  }

  // If we are getting the first child, we should start by filling the breeder with selected parents and breeding to
  // make ALL the children.
  if (child_index == 0)
  {
    fill();
    breed();
  }

  // Return the referred to child, updating the child_index to point to the next
  return std::move(child_[child_index++]);
}

template <class Solution, class Selector>
void OrgyBreeder<Solution, Selector>::reset()
{
  child_index = 0;
}

template <class Solution, class Selector>
void OrgyBreeder<Solution, Selector>::setSize(int newSize)
{
  size_ = newSize;
  child_ = std::vector<std::unique_ptr<Solution> >(size_);
  child_index = 0;
}

template <class Solution, class Selector>
void OrgyBreeder<Solution, Selector>::setMutationProb(double newProb)
{
  mutation_prob = newProb;
}

template <class Solution, class Selector>
void OrgyBreeder<Solution, Selector>::setCrossoverProb(double newProb)
{
  crossover_prob = newProb;
}

template <class Solution, class Selector>
void OrgyBreeder<Solution, Selector>::fill()
{
  // Creates copies of the selected parent solutions. This, followed by crossover, is somewhat inefficient.
  // Also inefficient to create and then copy.
  for (auto& child : child_)  // (child_ is the array of children, child is one of them. Rename?)
  {
    child = std::make_unique<Solution>(*selector_.select());  // Creates a new solution, using copy constructor.
  }
}

template <class Solution, class Selector>
void OrgyBreeder<Solution, Selector>::breed()
{
  // Shuffle the parents, to enable the simple application of crossover using selection without replacement
  std::shuffle(child_.begin(), child_.end(), utils::rand::rng);

  // Now we can breed in pairs, starting from the top.
  for (int i = 0; i < size_ - 1; i += 2)   // Must have -1 in case there is an odd number
  {
    if (utils::rand::rand01() < crossover_prob)
    {
      crossover(*child_[i], *child_[i + 1]);
    }
    child_[i]->mutate(mutation_prob);
    child_[i + 1]->mutate(mutation_prob);
  }

  if (size_ % 2 == 1)
  {
    // If there is an odd number, the last solution does not undergo crossover, but mutation may occur
    child_[size_ - 1]->mutate(mutation_prob);
  }
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution, class Selector>
GPBreeder<Solution, Selector>::GPBreeder(Selector& selector):
selector_(selector)  // (Mutation and crossover probabilities are uninitialized!)
{
}

template <class Solution, class Selector>
std::unique_ptr<Solution> GPBreeder<Solution, Selector>::child()
{
  // If we are getting the first child, first create the children required
  if (child_index == 0)
  {
    performing_crossover = (utils::rand::randDouble(0, 1) < crossover_prob);
    fill();
    breed();
  }

  // Return the referred to child, updating the child_index to point to the next
  auto oldIndex = child_index;
  if (performing_crossover && child_index == 0)
  {
    // Crossover has just been applied, so update child_index to refer to second child
    child_index = 1;
  }
  else
  {
    // Set child_index to zero to ensure that a new child is created next time
    child_index = 0;
  }

  return std::move(child_[oldIndex]);
}

template <class Solution, class Selector>
void GPBreeder<Solution, Selector>::reset()
{
  // Resets the breeder, between iterations, to ensure that new children are created the next time 'child' is called.
  child_index = 0;
}

template <class Solution, class Selector>
void GPBreeder<Solution, Selector>::setCrossoverProb(double newProb)
{
  crossover_prob = newProb;
}

template <class Solution, class Selector>
void GPBreeder<Solution, Selector>::fill()
{
  // Creates copies of the selected parent solutions. This, followed by crossover, is somewhat inefficient.
  // Also inefficient to create and then copy.
  child_[0] = std::make_unique<Solution>(*selector_.select());  // Creates a new solution, using copy constructor.
  if (performing_crossover)
  {
    child_[1] = std::make_unique<Solution>(*selector_.select());
  }
}

template <class Solution, class Selector>
void GPBreeder<Solution, Selector>::breed()
{
  if (performing_crossover)
  {
    crossover(*child_[0], *child_[1]);
  }
  else
  {
    child_[0]->mutate();
  }
}

#endif
