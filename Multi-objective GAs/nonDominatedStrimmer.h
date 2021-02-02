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

#ifndef NON_DOMINATED_STRIMMER
#define NON_DOMINATED_STRIMMER

#include <vector>

// Class to handle the distance matrices and the comparison of items based first on the distance of the nearest items,
// then the second nearest, and so on. Used when the number of non-dominated solutions exceeds the archive size in
// SPEA2.
class NonDominatedStrimmer
{
  friend class DistanceIterator;
public:
  explicit NonDominatedStrimmer(int numNonDominated);
  void addDistance(int i, int j, double dist);
  void initializeIndices();

  int mostCrowded() const;
  void remove(int toRemove);

private:
  bool moreCrowded(int lhs, int rhs) const;

  class DistanceIterator
  {
  private:
    const NonDominatedStrimmer* comp_;

    int base_;  // Remains constant on application of ++
    int i_;     // Moves to the next valid index on application of ++

  public:
    DistanceIterator(const NonDominatedStrimmer& comp, int base, int i);
    DistanceIterator(const DistanceIterator& rhs);

    bool operator==(const DistanceIterator& rhs) const;
    bool operator!=(const DistanceIterator& rhs) const;

    int operator*() const;
    DistanceIterator& operator++();
    DistanceIterator  operator++(int);
  };

  DistanceIterator begin(int base) const;
  DistanceIterator end(int base) const;


private:
  // Distance matrix
  int num_items;
  std::vector<std::vector<double> > distances_;

  // Each index vector holds the indices of the other objects, sorted on distance.
  // However, this sorting takes place in a lazy manner
  std::vector<std::vector<int> > index_by_distance;

  // Vector of bools, with false indicating that the object has been removed.
  std::vector<bool> available_;
};

#endif
