// MOQCD1 (Multi-objective quantum circuit discovery - algorithm 1.)
//
// Copyright (c) 2020  Alan Reynolds (alanreynolds9@gmail.com)
//
// This file is part of MOQCD1.
//
// MOQCD1 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
// MOQCD1 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with MOQCD1.  If not, see
// <http://www.gnu.org/licenses/>.

//----------------------------------------------------------------------------------------------------------------------

#ifndef CACHE_H
#define CACHE_H

#include <vector>
#include <memory>
#include <unordered_map>
#include <shared_mutex>
#include "circuit.h"

// First a couple of classes providing function objects for the hash function and equivalence function used in the
// operation of the cache.

class CircuitHash
{
public:
  size_t operator()(const circuit::Circuit& circuit) const;
};


class CircuitsEquivalent
{
public:
  bool operator()(const circuit::Circuit& lhs, const circuit::Circuit& rhs) const;
};



// Next, a class that holds the statistics regarding circuit quality, numerical optimization results and history, etc.,
// that is cached for each simplified (canonical) circuit structure visited, along with the simplified (canonical)
// circuit itself.

class CachedStats
{
public:
  const circuit::Circuit& circuit() const;
  const std::vector<double>& error() const;
  double primaryError() const;
  std::shared_ptr<CachedStats> reduced() const;  // Follow this for statistics for the reduced circuit.

  void setCircuit(const circuit::Circuit& circuit);
  void setEvaluated(const std::vector<double>& error);
  void setParameters(const std::vector<double>& values);
  void setReduced(std::shared_ptr<CachedStats> reducedStats);

  // While it's tempting just to make the mutex public, we instead provide the following two functions to control thread
  // access to the CachedStats object.
  std::shared_lock<std::shared_timed_mutex> getReadersLock() const;
  std::unique_lock<std::shared_timed_mutex> getWritersLock() const;

private:
  std::unique_ptr<circuit::Circuit> circuit_;

  // The errors associated with the 'best' circuit with the supplied circuit structure, i.e. the result of the numerical
  // optimization.
  std::vector<double> error_;

  // Pointer to the statistics for the 'reduced' version of the circuit, i.e. the version where gates with angles near
  // zero removed and those with angles near 'special' values replaced with cheaper alternatives, if reduction is
  // possible. Otherwise just the nullptr. (Note that the angle values stored in circuit_ should match the reduction
  // performed. Reapplying numerical optimization should result in this pointer first being reset to nullptr, then set
  // to point to the reduction of the newly optimized circuit, if a reduction exists.)
  std::shared_ptr<CachedStats> reduced_;
  mutable std::shared_timed_mutex mtx_;  // All I want is a shared_mutex, but the Mac's compiler can't find it!
};



// Finally, the cache itself.

class CircuitCache
{
public:
  CircuitCache();

  // Function that finds a circuit in the cache, if present. If not, it returns a nullptr.
  std::shared_ptr<CachedStats> find(const circuit::Circuit& circuit) const;

  // Two functions that add links from circuits to cached statistics. The first is used when the simplified circuit is
  // already in the cache, but new circuits are discovered that simplify to it. The new circuits are linked to the
  // provided statistics object. The second completes a new statistics object by inserted a copy of the new simplified
  // circuit. In both cases, we ensure that pointers added to the cache are not null and do point to something.
  int size() const;
  void add(const std::vector<circuit::Circuit>& circuits, std::shared_ptr<CachedStats> stats);
  std::shared_ptr<CachedStats> add(const std::vector<circuit::Circuit>& circuits, const circuit::Circuit& simplified,
                                   std::shared_ptr<CachedStats> emptyStats);
  int numHits() const;

private:
  std::shared_ptr<CachedStats> find_unthreaded(const circuit::Circuit& circuit) const;
  void add_unthreaded(const std::vector<circuit::Circuit>& circuits, std::shared_ptr<CachedStats> stats);

private:
  using HashTable = std::unordered_map<circuit::Circuit, std::shared_ptr<CachedStats>, CircuitHash, CircuitsEquivalent>;

  HashTable cache_;
  mutable std::shared_timed_mutex mtx_;  // All I want is a shared_mutex, but the Mac's compiler can't find it!
  mutable std::atomic<int> num_hits = 0;
};


#endif  // CACHE_H
