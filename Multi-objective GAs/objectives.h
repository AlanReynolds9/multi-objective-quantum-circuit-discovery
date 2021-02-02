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

#ifndef OBJECTIVES_H
#define OBJECTIVES_H

#include <cassert>
#include <vector>

// Objectives class

class Objectives
{
  // (We need to think about the interface of this class. Functions such as the new 'scaleObjective' may be preferable
  // to the lower level functions such as 'lowerBound'. We also need to have a sensible way of dealing with unknown
  // objective bounds.)
public:
  template <class Problem> explicit Objectives(const Problem& prob);

  int numObj() const;
  bool maximized(int objNum) const;
  double lowerBound(int objNum) const;  // Upper and lower bounds have only been used to calculate objective ranges...
  double upperBound(int objNum) const;  // ...thus far. Also, such ranges may be unknown for some problems.
  double scaleObjective(int objNum, double value) const;

private:
  int num_obj;
  std::vector<bool> maximized_;     // Currently need only the range of the objectives. However, we will keep the...
  std::vector<double> lower_bound;  // ...rest of this information in case, for example, we decide to implement...
  std::vector<double> upper_bound;  // ...simple roulette selection on one objective.
};


// Comparator classes for comparing solutions by objective

template <class Solution, template <typename> class Compare>
class CompareByObjective
{
  // Function object that determines whether one solution has a smaller objective value than another for a selected
  // objective. Used by NSGA II crowding and by the Store classes. As with the rest of the code, assumes that double is
  // used as the objective type.
public:
  explicit CompareByObjective(int objNum);
  void selectObjective(int objNum);
  bool operator()(const Solution& lhs, const Solution& rhs) const;

  template <class Pointer>
  bool operator()(Pointer lhs, Pointer rhs) const;

private:
  int obj_num;
};

template <class Solution>
class LexBetterByObjective
{
  // Same as CompareByObjective, but, in the case of a tie, proceeds to compare against the next objective. Used in some
  // of the NSGA II evaluator code. Assumes that double is the objective type. (Also assumes that a 'tie' according to
  // the comparator implies equal objective values.)
public:
  LexBetterByObjective(const Objectives& objectives, int objNum);
  bool operator()(const Solution& lhs, const Solution& rhs) const;

  template <class Pointer>
  bool operator()(Pointer lhs, Pointer rhs) const;

private:
  const Objectives* objectives_;
  int obj_num;
};

//----------------------------------------------------------------------------------------------------------------------

template <class Problem>
Objectives::Objectives(const Problem& prob) :
num_obj(prob.numObj()),
maximized_(prob.numObj()),
lower_bound(prob.numObj()),
upper_bound(prob.numObj())
{
  for (int i = 0; i < numObj(); ++i)
  {
    maximized_[i] = prob.objMaximized(i);
    lower_bound[i] = prob.objLowerBound(i);
    upper_bound[i] = prob.objUpperBound(i);
  }
}

inline int Objectives::numObj() const
{
  return num_obj;
}

inline bool Objectives::maximized(int objNum) const
{
  return maximized_[objNum];
}

inline double Objectives::lowerBound(int objNum) const
{
  return lower_bound[objNum];
}

inline double Objectives::upperBound(int objNum) const
{
  return upper_bound[objNum];
}

inline double Objectives::scaleObjective(int objNum, double value) const
{
  // Takes a raw objective value and converts to an objective ranging from zero to one, to be maximized.
  // Warning: Should only be used when true upper and lower bounds for the objective values are known, otherwise values
  // outside [0,1] may be returned.
  double scaledValue;
  if (maximized_[objNum])
  {
    scaledValue = value - lower_bound[objNum];
  }
  else
  {
    scaledValue = upper_bound[objNum] - value;
  }

  scaledValue /= upper_bound[objNum] - lower_bound[objNum];

  return scaledValue;
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution, template <typename> class Compare>
CompareByObjective<Solution, Compare>::CompareByObjective(int objNum) :
obj_num(objNum)
{
}

template <class Solution, template <typename> class Compare>
void CompareByObjective<Solution, Compare>::selectObjective(int objNum)
{
  obj_num = objNum;
}

template <class Solution, template <typename> class Compare>
bool CompareByObjective<Solution, Compare>::operator()(const Solution& lhs, const Solution& rhs) const
{
  return Compare<double>()(lhs.objective(obj_num), rhs.objective(obj_num));
}

template <class Solution, template <typename> class Compare>
template <class Pointer>
bool CompareByObjective<Solution, Compare>::operator()(Pointer lhs, Pointer rhs) const
{
  return (*this)(*lhs, *rhs);
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution>
LexBetterByObjective<Solution>::LexBetterByObjective(const Objectives& objectives, int objNum) :
objectives_(&objectives),
obj_num(objNum)
{
  assert(obj_num >= 0 && obj_num < static_cast<int>(objectives_->numObj()));  // (Would like to lose the cast.)
}

template <class Solution>
bool LexBetterByObjective<Solution>::operator()(const Solution& lhs, const Solution& rhs) const
{
  int objToCompare = obj_num;
  do
  {
    if (lhs.objective(objToCompare) != rhs.objective(objToCompare))
    {
      if (objectives_->maximized(objToCompare))
      {
        return lhs.objective(objToCompare) > rhs.objective(objToCompare);
      }
      else
      {
        return lhs.objective(objToCompare) < rhs.objective(objToCompare);
      }
    }

    ++objToCompare;
    if (objToCompare == objectives_->numObj())
    {
      objToCompare = 0;
    }
  }
  while (objToCompare != obj_num);

  return false;
}

template <class Solution>
template <class Pointer>
bool LexBetterByObjective<Solution>::operator()(Pointer lhs, Pointer rhs) const
{
  return (*this)(*lhs, *rhs);
}


#endif
