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

#include "cache.h"
#include "circuit.h"

using std::shared_ptr;
using std::make_unique;
using std::make_shared;
using std::vector;
using circuit::Circuit;


size_t CircuitHash::operator()(const Circuit& circuit) const
{
  return circuit.hash();
}

//----------------------------------------------------------------------------------------------------------------------

bool CircuitsEquivalent::operator()(const Circuit& lhs, const Circuit& rhs) const
{
  return lhs.equivalentStructure(rhs);
}


//-------------------------------------------------------------------------------

const vector<double>& CachedStats::error() const
{
  return error_;
}


double CachedStats::primaryError() const
{
  return error_[0];
}


shared_ptr<CachedStats> CachedStats::reduced() const
{
  // Pointer to the CachedStats of the reduced circuit, i.e. that circuit that results from eliminating gates with angle
  // parameters close to zero and reducing those with special (or near to special) angles to cheaper gates.
  // Re-optimization can potentially lead to a chain of such links.
  return reduced_;
}


void CachedStats::setCircuit(const circuit::Circuit& circuit)
{
  circuit_ = make_unique<circuit::Circuit>(circuit);
}


void CachedStats::setEvaluated(const std::vector<double>& error)
{
  error_ = error;
}


void CachedStats::setParameters(const vector<double>& values)
{
  circuit_->setParameters(values);
}


void CachedStats::setReduced(shared_ptr<CachedStats> reducedStats)
{
  reduced_ = reducedStats;
}


const Circuit& CachedStats::circuit() const
{
  return *circuit_;
}


std::shared_lock<std::shared_timed_mutex> CachedStats::getReadersLock() const
{
  return std::shared_lock<std::shared_timed_mutex>(mtx_);
}


std::unique_lock<std::shared_timed_mutex> CachedStats::getWritersLock() const
{
  return std::unique_lock<std::shared_timed_mutex>(mtx_);
}

//-------------------------------------------------------------------------------

CircuitCache::CircuitCache()
{
  cache_.reserve(1000000);  // MAGIC NUMBER: I don't know how many circuits will be generated, but it will be a lot!
}


int CircuitCache::size() const
{
  return static_cast<int>(cache_.size());
}


shared_ptr<CachedStats> CircuitCache::find_unthreaded(const circuit::Circuit& circuit) const
{
  // Find the CachedStats associated with a given circuit in the cache. If the circuit is not in the cache, return the
  // null_ptr. This version is not thread-safe and is therefore private, but it is called by the thread-safe version and
  // by one of the cacheing functions.
  auto i = cache_.find(circuit);
  if (i == cache_.end())
  {
    return nullptr;
  }
  ++num_hits;
  return i->second;
}


shared_ptr<CachedStats> CircuitCache::find(const Circuit& circuit) const
{
  std::shared_lock<std::shared_timed_mutex> lck{mtx_};  // Get lock for reading.
  return find_unthreaded(circuit);
}


void CircuitCache::add_unthreaded(const std::vector<circuit::Circuit> &circuits, std::shared_ptr<CachedStats> stats)
{
  // Cache links from new (unsimplified) circuits to an existing CachedStats object for some previously seen simplified
  // circuit. This version is not thread-safe and is therefore private, but it is called by thread-safe versions.
  for (const auto& circuit : circuits)
  {
    // Note that, on occasions another thread might insert the same circuit into the cache first. When this happens, the
    // other thread should associate the circuit with the same 'stats' pointer, so when the line below does nothing,
    // that is precisely what we want to happen. (I've done some limited testing of this.)
    cache_.insert(make_pair(circuit, stats));
  }
}


void CircuitCache::add(const vector<Circuit>& circuits, shared_ptr<CachedStats> stats)
{
  // Cache links from new (unsimplified) circuits to an existing CachedStats object for the previously seen simplified
  // circuit. Thread-safe.
  std::unique_lock<std::shared_timed_mutex> lck{mtx_};  // Get lock for writing.
  add_unthreaded(circuits, stats);
}


shared_ptr<CachedStats> CircuitCache::add(const vector<Circuit>& circuits, const Circuit& simplified,
                                          shared_ptr<CachedStats> emptyStats)
{
  // This functions usually puts the new simplified circuit into the emptyStats object and then caches links from all
  // the circuits in 'circuits', i.e. the simplified circuit and all the unsimplified circuits that simplified to it, to
  // the provided CachedStats object. However, there are situations where a CachedStats object for the circuit is
  // already present in the cache, in which case we cache links to that instead.
  std::unique_lock<std::shared_timed_mutex> lck{mtx_};  // Get lock for writing
  auto stats = find_unthreaded(simplified);
  if (stats)
  {
    // The code currently gets here in two situations.
    // 1. The circuit has been simplified using Circuit::simplifyCircuit(), which ignores existing elements in the cache
    //    and calls this add() function, regardless of whether the simplified circuit has been seen before.
    // 2. Two threads discover the same simplified circuit simultaneously. One successfully completes and adds
    //    emptyStats using the code at the bottom of this function. The other ends here. (Rare, though tests reveal that
    //    this occurs in perhaps 10% of runs of 300 generations with population of 100.)
    add_unthreaded(circuits, stats);
    return stats;
  }

  // The typical situation where we are adding a new simplified circuit to the cache. Insert the new simplified circuit
  // into the empty CachedStats object and then cache links to it from each of the circuits.
  emptyStats->setCircuit(simplified);
  add_unthreaded(circuits, emptyStats);
  return emptyStats;
}


int CircuitCache::numHits() const
{
  return num_hits;
}
