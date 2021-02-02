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

// Function objects for calculating distances based on objective values

#ifndef DISTANCES_H
#define DISTANCES_H

#include <limits>
#include <exception>
#include "objectives.h"

// (For now this is ok. However, objectives usually need to be scaled in order to have sensible distances. Such scaling
// might be static, if objective ranges are provided or dynamic if scaling is perfomed with respect to the objective
// values present in the population.)

// First a templated function to find the distance between two solutions, by objective.
template <class Solution, int d>
class ObjectiveNorm
{
public:
  explicit ObjectiveNorm(const Objectives& objectives) :
  objectives_(objectives)
  {
    // (Would prefer some way of using a compile time assertion. Alternatively, some way to automatically use an
    // alternative.)
    for (auto obj = 0; obj < objectives_.numObj(); ++obj)
    {
      if (objectives_.upperBound(obj) - objectives_.lowerBound(obj) == 0)
      {
        throw std::domain_error("Attempting to use the ObjectiveNorm class when the objective ranges are unknown.");
      }
    }
  }

  double operator()(const Solution& lhs, const Solution& rhs) const
  {
    // Calculate the distance according to the d-norm, i.e. dth root of (x1-x2)^d + (y1-y2)^d
    double distance = 0;
    for (auto objNum = 0; objNum < objectives_.numObj(); ++objNum)
    {
      // Find the difference between the objective values for this objective and
      // normalise according to the min and max values supplied by the Solution class
      double difference = abs(lhs.objective(objNum) - rhs.objective(objNum));
      difference /= objectives_.upperBound(objNum) - objectives_.lowerBound(objNum);

      // Adjust distance according to the d'th power of the objective difference
      distance += pow(difference, d);
    }

    return pow(distance, 1.0/d);
  }

private:
  const Objectives& objectives_;
};

template <class Solution>
class ObjectiveNorm<Solution, 1>
{
  // Specialization for the Manhattan metric, for efficiency reasons
public:
  explicit ObjectiveNorm(const Objectives& objectives) :
  objectives_(objectives)
  {
    // (Would prefer some way of using a compile time assertion. Alternatively, some way to automatically use an
    // alternative.)
    for (auto obj = 0; obj < objectives_.numObj(); ++obj)
    {
      if (objectives_.upperBound(obj) - objectives_.lowerBound(obj) == 0)
      {
        throw std::domain_error("Attempting to use the ObjectiveNorm class when the objective ranges are unknown.");
      }
    }
  }

  double operator()(const Solution& lhs, const Solution& rhs) const
  {
    double distance = 0;
    for (auto objNum = 0; objNum < objectives_.numObj(); ++objNum)
    {
      // Find the difference between the objective values for this objective and
      // normalise according to the min and max values supplied by the Solution class
      double difference = abs(lhs.objective(objNum) - rhs.objective(objNum));
      difference /= objectives_.upperBound(objNum) - objectives_.lowerBound(objNum);

      // Adjust distance by adding the objective difference
      distance += difference;
    }

    return distance;
  }

private:
  const Objectives& objectives_;
};

template <class Solution>
class ObjectiveNorm<Solution, 2>
{
  // Specialization for the Euclidean metric, for efficiency reasons only
public:
  explicit ObjectiveNorm(const Objectives& objectives) :
  objectives_(objectives)
  {
    // (Would prefer some way of using a compile time assertion. Alternatively, some way to automatically use an
    // alternative.)
    for (auto obj = 0; obj < objectives_.numObj(); ++obj)
    {
      if (objectives_.upperBound(obj) - objectives_.lowerBound(obj) == 0)
      {
        throw std::domain_error("Attempting to use the ObjectiveNorm class when the objective ranges are unknown.");
      }
    }
  }

  double operator()(const Solution& lhs, const Solution& rhs) const
  {
    double distance = 0;
    for (auto objNum = 0; objNum < objectives_.numObj(); ++objNum)
    {
      // Find the difference between the objective values for this objective and
      // normalise according to the min and max values supplied by the Solution class
      double difference = lhs.objective(objNum) - rhs.objective(objNum);
      difference /= objectives_.upperBound(objNum) - objectives_.lowerBound(objNum);

      // Adjust distance by adding the square of the objective difference
      distance += difference * difference;
    }

    return sqrt(distance);
  }

private:
  const Objectives& objectives_;
};

#endif
