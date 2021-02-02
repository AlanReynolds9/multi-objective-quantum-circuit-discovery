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

#include <queue>
#include <unordered_set>
#include <LBFGS.h>
#include "cache.h"
#include "circuit.h"
#include "replacements.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
using std::make_pair;
using std::unique_ptr;
using std::shared_ptr;
using std::make_shared;
using std::make_move_iterator;
using std::vector;
using utils::rand::randInt;
using utils::rand::rand01;
using utils::rand::randGeometric;

using Eigen::VectorXd;
using Eigen::Map;
using namespace LBFGSpp;


void deletePointers(vector<const Gate*>& pointers)
{
  for (auto pointer : pointers)
  {
    delete pointer;
  }
}


// Vector conversion functions and the error and gradient function for LBFGS++.
VectorXd vectorToVectorXd(vector<double>& in)
{
  // Warning: This merely LINKS to the data in the std::vector. As a result, the std::vector will be changed whenever
  // the VectorXd is. (Hence the input parameter is not const.
  return Map<VectorXd>(in.data(), in.size());
}


vector<double> vectorXdToVector(const VectorXd& in)
{
  return vector<double>(in.data(), in.data() + in.size());
}


// GateCost class

GateCost::GateCost()
{
}


GateCost::GateCost(const vector<long>& coeff) :
coeff_(coeff)
{
}


GateCost::GateCost(vector<long>&& coeff) :
coeff_(std::move(coeff))
{
}


GateCost::GateCost(std::initializer_list<long> coeff) :
coeff_(coeff)
{
}


long GateCost::cost(int numControls) const
{
  long power{1};
  long total{0};
  for (auto c : coeff_)
  {
    total += c * power;
    power *= numControls;
  }

  return total;
}


void GateCost::output(std::ostream& out) const
{
  if (coeff_.empty())
  {
    out << "uninitialized" << endl;
    return;
  }
  out << coeff_[0];
  for (auto i = 1; i < coeff_.size(); ++i)
  {
    out << " + " << coeff_[i] << "c";
    if (i > 1)
    {
      out << "^i";
    }
  }
}


ostream& operator<<(ostream& out, const GateCost& gateCost)
{
  gateCost.output(out);
  return out;
}


namespace circuit
{
  class OverallErrorAndGradient
  {
    // Class for use by the numerical optimizer LBFGS++. (Hence the odd datatypes.)
  public:
    OverallErrorAndGradient(Solution& solution) :
    solution_(solution)
    {
    }

    double operator()(const VectorXd& parameters, VectorXd &grad)
    {
      // Copy the input parameters into a std::vector.
      vector<double> parameterVector = vectorXdToVector(parameters);

      // Copy the parameters into the circuit.
      solution_.setParameters(parameterVector);

      // Find the overall error and gradient and put the gradient into gradOut.
      solution_.evaluate_error_and_gradient();
      vector<double> gradient = solution_.gradient_;  // vectorToVectorXd needs a non-const (hence non-temporary)...
      grad = vectorToVectorXd(gradient);              // ...std::vector to work with.

      // Return the overall error.
      evals++;
      return solution_.primary_error();
    }

  private:
    Solution& solution_;

  public:
    static int evals;
  };

  int OverallErrorAndGradient::evals = 0;

  void outputEvals(ostream& out)
  {
    out << "Evals = " << OverallErrorAndGradient::evals << endl;
  }

  //----------------------------------------------------------------------------------------------------------------

  CircuitContext::CircuitContext(const vector<GateCreator>& gateCreator,
                                 const vector<PermittedControls>& permittedControls,
                                 const vector<GateCost>& gateCosts, int numQbits, double meanRandomLength) :
  num_permitted_gate_types(static_cast<int>(gateCreator.size())),
  num_permitted_gate_options(0),  // Incremented in the function body.
  gate_creator(gateCreator),
  gate_type_available(constants::totalNumGateTypes, vector<bool>(numQbits, false)),
  gate_cost(constants::totalNumGateTypes),
  gate_option_base_id(constants::totalNumGateTypes, 0),
  num_qbits(numQbits),
  qbit_input_options(numQbits, QbitOptions::varies),
  mean_random_length(meanRandomLength)
  {
    assert(permittedControls.size() == num_permitted_gate_types);
    assert(gateCosts.size() == num_permitted_gate_types);

    for (auto i = 0; i < num_permitted_gate_types; ++i)
    {
      // Get a gate of the permitted type and the gate type id.
      auto gate = gate_creator[i](*this);  // UNPLEASANT: Gate created (with broken option ID) before CircuitContext
      int id = gate->gateTypeId();         //             is completely constructed.

      // Set gate costs
      gate_cost[id] = gateCosts[i];

      // Set gate availability. Index of gate_type_available[] is 0, 1 or 2 but means, none, one and many respectively.
      // Gates that don't have controls will look only at gate_type_available[0].
      switch (permittedControls[i])
      {
        case PermittedControls::none:
          gate_type_available[id][0] = true;
          break;

        case PermittedControls::one:
          gate_type_available[id][1] = true;
          break;

        case PermittedControls::notMany:
          gate_type_available[id][0] = gate_type_available[id][1] = true;
          break;

        case PermittedControls::many:
          make_available_with_many_controls(id);
          break;

        case PermittedControls::notOne:
          gate_type_available[id][0] = true;
          make_available_with_many_controls(id);
          break;

        case PermittedControls::atLeastOne:
          gate_type_available[id][1] = true;
          make_available_with_many_controls(id);
         break;

        case PermittedControls::any:
          gate_type_available[id][0] = gate_type_available[id][1] = true;
          make_available_with_many_controls(id);
          break;

        default:
          throw std::domain_error("Parameter 'permittedControls' in CircuitContext::CircuitContext contains an invalid "
                                  "element.");
      }
    }

    // Remove availability of gate types that are redundant, i.e. that can be performed by a cheaper gate, or one of the
    // same cost that is more flexible.
    remove_redundancies();

    // Determine the base 'option ID' for each available gate type. Each combination of gate type, target and controls
    // gets a unique 'option ID', used for hashing. (Must occur AFTER setting gate_type_available the the gate types.)
    for (auto i = 0; i < num_permitted_gate_types; ++i)
    {
      // Get a gate of the permitted type and the gate type id.
      // (Use gate_creator, NOT gateCreator, which will not have had redundant gates removed.)
      auto gate = gate_creator[i](*this);  // UNPLEASANT: Gate created (with broken option ID) before CircuitContext
      int id = gate->gateTypeId();         //             is completely constructed.

      // Calculate base option ID.
      gate_option_base_id[id] = num_permitted_gate_options;
      num_permitted_gate_options += gate->numQbitOptions();
    }
  }


  int CircuitContext::numPermittedGateOptions() const
  {
    return num_permitted_gate_options;
  }


  bool CircuitContext::gateTypeAvailable(int gateTypeId, int numControls) const
  {
    // For use when we want to avoid creating a gate just to check availability. The Gate classes contain static
    // 'available' routines that call this function, with uncontrolled gate classes such as SwapGate calling it with
    // numControls set to zero.
    assert(numControls >= 0);
    assert(numControls < num_qbits);
    return gate_type_available[gateTypeId][numControls];
  }


  bool CircuitContext::gateTypeAvailable(const Gate& gate, int numControls) const
  {
    // The Gate classes contain their own 'available()' routines that call this function. Uncontrolled gate classes,
    // such as SwapGate, just call this function with numControls set to zero.
    return gate_type_available[gate.gateTypeId()][numControls];
  }


  long CircuitContext::gateCost(int gateTypeId, int numControls) const
  {
    // For use when we want to avoid creating a gate just to check availability.
    return gate_cost[gateTypeId].cost(numControls);
  }


  long CircuitContext::gateCost(const Gate& gate, int numControls) const
  {
    return gate_cost[gate.gateTypeId()].cost(numControls);
  }


  size_t CircuitContext::gateOptionBaseId(const Gate& gate) const
  {
    return gate_option_base_id[gate.gateTypeId()];
  }


  int CircuitContext::numQbits() const
  {
    return num_qbits;
  }


  void CircuitContext::setQbitInputOptions(int qbit, QbitOptions options)
  {
    qbit_input_options[qbit] = options;
  }


  QbitOptions CircuitContext::qbitInputOptions(int qbit) const
  {
    return qbit_input_options[qbit];
  }


  double CircuitContext::meanRandomLength() const
  {
    return mean_random_length;
  }


  unique_ptr<Gate> CircuitContext::randomGate() const
  {
    // Creates a random Gate, permitted by the CircuitContext, i.e. of a type in the set selected by the user.

    // Get a gate of random type.
    auto gateType = randInt(0, num_permitted_gate_types);
    auto gate = create_selected_gate(gateType);

    // Make it a random gate of that type.
    gate->random();

    // Compiler should use 'return value optimization' or 'copy elision'. No need for std::move.
    return gate;
  }


  void CircuitContext::output(ostream& out) const
  {
    out << "Number of qbits = " << numQbits() << endl;
    out << "Mean length of a random solution = " << meanRandomLength() << endl;
    out << "Number of available gate types = " << num_permitted_gate_types << endl;
    out << "Number of available gate options (type/target/controls) = " << num_permitted_gate_options << endl;
    out << "Gate types available:" << endl;
    for (auto i = 0; i < num_permitted_gate_types; ++i)
    {
      // Get the gate type id and name.
      auto gate = gate_creator[i](*this);
      int id = gate->gateTypeId();
      out << "  " << gate->name();

      // Output which types of controls are available.
      bool allMany = true;
      bool anyMany = false;
      for (int numControls = 2; numControls < numQbits(); ++numControls)
      {
        if (gate_type_available[id][numControls])
        {
          anyMany = true;
        }
        else
        {
          allMany = false;
        }
      }
      if (!anyMany)
      {
        if (gate_type_available[id][0])
        {
          if (gate_type_available[id][1])
          {
            out << " with zero or one control.";
          }
          else
          {
            out << " with no controls.";
          }
        }
        else
        {
          if (gate_type_available[id][1])
          {
            out << " with one control.";
          }
          else
          {
            out << ", but there is no permitted option for controls!";
          }
        }
      }
      else if (allMany)
      {
        if (gate_type_available[id][0])
        {
          if (gate_type_available[id][1])
          {
            out << " with any number of controls.";
          }
          else
          {
            out << " with no controls, or at least two!";
          }
        }
        else
        {
          if (gate_type_available[id][1])
          {
            out << " with at least one control.";
          }
          else
          {
            out << " with at least two controls!";
          }
        }
      }
      else
      {
        out << " with number of controls selected from the set \{ ";
        bool first = true;
        for (int numControls = 0; numControls < numQbits(); ++numControls)
        {
          if (gate_type_available[id][numControls])
          {
            if (first)
            {
              out << numControls;
              first = false;
            }
            else
            {
              out << ", " << numControls;
            }
          }
        }
      }

      out << " Cost = " << gate_cost[id] << "." << " Gate option IDs start from " << gate_option_base_id[id] << "."
          << endl;
    }
  }


  void CircuitContext::make_available_with_many_controls(int gateTypeId)
  {
    for (int numControls = 2; numControls < num_qbits; ++numControls)
    {
      gate_type_available[gateTypeId][numControls] = true;
    }
  }


  void CircuitContext::remove_redundancies()
  {
    // Since specific phase shift gates such as PiByEight, PhaseInv and Z can be implemented as an ArbitraryPhase, we
    // remove availability of these gates when the ArbitraryPhase is not more expensive.
    // Note that we could write this function without requiring the creation of the various gates by hard coding the
    // gateTypeId's. However, this could potentially break if we decide to add new gates or reorder them.

    // Find the gate type ID for each of the involved gates.
    int arbPhaseId = ArbitraryPhase(*this).gateTypeId();
    int piByEightId = PiByEight(*this).gateTypeId();
    int piByEightInvId = PiByEightInv(*this).gateTypeId();
    int phaseId = PhaseGate(*this).gateTypeId();
    int phaseInvId = PhaseInv(*this).gateTypeId();
    int zId = ZGate(*this).gateTypeId();
    int zRotId = ZRotation(*this).gateTypeId();

    // Eliminate redundant options amongst the phase type gates, e.g. PhaseGate with no controls when ArbitraryPhase
    // with no controls is available at the same cost.
    vector<std::pair<int, string> > phaseTypes{{piByEightId, "PiByEight"}, {piByEightInvId, "PiByEightInv"},
                                               {phaseId, "PhaseGate"}, {phaseInvId, "PhaseInv"}, {zId, "ZGate"}};
    for (auto numControls = 0; numControls < numQbits(); ++numControls)
    {
      if (gate_type_available[arbPhaseId][numControls])
      {
        for (auto [id, name] : phaseTypes)
        {
          if (gate_type_available[id][numControls] &&
              gate_cost[arbPhaseId].cost(numControls) <= gate_cost[id].cost(numControls))
          {
            cout << "Warning: Both ArbitraryPhase and " << name << " are available with " << numControls << "control";
            if (numControls != 1)
            {
              cout << "s";
            }
            cout << ", with ArbitraryPhase no more expensive than " << name << ". ";
            cout << "Removing availability of " << name << " with " << numControls << "control";
            if (numControls != 1)
            {
              cout << "s";
            }
            cout << "." << endl;

            gate_type_available[id][numControls] = false;
          }
        }
      }
    }

    // Eliminate redundant uncontrolled ZRotation.
    // NOTE: While we used to consider a controlled ZRotation to be redundant if it could be replaced, at reduced cost,
    // by two ArbitraryPhase gates, this is more complicated to detect now that gate cost depends on the number of
    // controls. Moreover, it is unclear whether it is a good idea to mark such gates as redundant anyway, since
    // circuits with fewer gates are easier to find. A better way of handling this might be to add single ZRotation
    // simplification to the code?
    if (gate_type_available[zRotId][0])
    {
      if (gate_type_available[arbPhaseId][0] && gate_cost[arbPhaseId].cost(0) <= gate_cost[zRotId].cost(0))
      {
        cout << "Warning: Both ArbitraryPhase and ZRotation are available without controls, with ArbitraryPhase no "
                "more expensive than ZRotation. Removing availability of ZRotation with zero controls." << endl;
        gate_type_available[zRotId][0] = false;
      }
    }

    // Eliminate gate types that no longer have any options available.
    int i{0};
    while (i < num_permitted_gate_types)
    {
      auto gate = gate_creator[i](*this);  // UNPLEASANT: Gate created (with broken option ID) before CircuitContext
      if (gate->numQbitOptions() == 0)     //             is completely constructed.
      {
        cout << "Warning: Gate type \"" << gate->name() << "\" has been removed entirely." << endl;
        gate_creator.erase(gate_creator.begin() + i);
        --num_permitted_gate_types;
      }
      else
      {
        ++i;
      }
    }
  }


  unique_ptr<Gate> CircuitContext::create_selected_gate(int gateTypeNum) const
  {
    return gate_creator[gateTypeNum](*this);
  }


  ostream& operator<<(ostream& out, const CircuitContext& circuitContext)
  {
    circuitContext.output(out);
    return out;
  }

  //------------------------------------------------------------------------------------------------------------------

  Problem::Problem(const vector<GateCreator>& gateCreator, const vector<PermittedControls>& permittedControls,
                   const vector<GateCost>& gateCosts, int numQbits, double meanRandomLength, double gateCostRange,
                   int maxLength, double targetError, long expectedReduction, bool randomizeParameters,
                   int iterationQuota, bool useCache) :
  circuit_context(gateCreator, permittedControls, gateCosts, numQbits, meanRandomLength),
  gate_cost_range(gateCostRange),
  max_length(maxLength),
  target_error(targetError),
  expected_reduction(expectedReduction),
  randomize_parameters(randomizeParameters),
  iteration_quota(iterationQuota),
  cache_(useCache ? std::make_shared<CircuitCache>() : nullptr)
  {
  }


  const CircuitContext& Problem::circuitContext() const
  {
    return circuit_context;
  }


  int Problem::numQbits() const
  {
    return circuit_context.numQbits();
  }


  void Problem::transferQbitInputOptionsToContext()
  {
    for (int qbit = 0; qbit < numQbits(); ++qbit)
    {
      circuit_context.setQbitInputOptions(qbit, qbitInputOptions(qbit));
    }
  }


  int Problem::maxLength() const
  {
    return max_length;
  }


  double Problem::targetError() const
  {
    return target_error;
  }


  long Problem::expectedReduction() const
  {
    return expected_reduction;
  }


  bool Problem::randomizeParameters() const
  {
    return randomize_parameters;
  }


  int Problem::iterationQuota() const
  {
    return iteration_quota;
  }


  shared_ptr<CircuitCache> Problem::cache() const
  {
    return cache_;
  }


  Front& Problem::front() const
  {
    return front_;
  }

void Problem::output(ostream& out) const
  {
    out << circuit_context;
    out << "Gate cost range = " << gate_cost_range << endl;
  }


  ostream& operator<<(ostream& out, const Problem& problem)
  {
    problem.output(out);
    return out;
  }

  //------------------------------------------------------------------------------------------------------------------

  Circuit::Circuit(const CircuitContext& context) :
  context_(&context),
  marked_state(-1)
  {
  }


  Circuit::Circuit(const CircuitContext& context, GateSequence&& gates) :
  context_(&context),
  gates_(std::move(gates)),
  marked_state(-2)
  {
    imprint_gates();  // Ensure gates get access to any data they need from their 'parent' circuit.
  }


  Circuit::Circuit(const Circuit& rhs) :
  context_(rhs.context_),
  marked_state(-3)  // Choosing an illegal value, rather than rhs.marked_state, catches bugs where a gate is incorrectly
  {                 // 'imprinted'.
    // Perform a deep copy of the gates_ vector. (We can't have the new Circuit having pointers into the old one!)

    gates_.reserve(rhs.length());
    for (const auto& gate : rhs.gates_)
    {
      gates_.push_back(gate->clone());
    }

    imprint_gates();  // Ensure gates get data from their new 'parent' circuit and not the old one.
  }


  Circuit::Circuit(Circuit&& rhs) :
  context_(std::exchange(rhs.context_, nullptr)),  // Don't have to exchange with nullptr - the (default) destructor...
  gates_(std::move(rhs.gates_)),                   // ...won't delete it anyway.
  marked_state(-4)  // Choosing an illegal value, rather than rhs.marked_state, catches bugs where a gate is...
  {                 // ...incorrectly imprinted.
    imprint_gates();  // Ensure gates get data they need from their new 'parent' circuit and not the old one.
  }


  Circuit& Circuit::operator=(const Circuit& rhs)
  {
    // Perform a deep copy of the gates_ vector. (We can't have two Circuits having pointers to the same gates!)

    // Check for self assignment.
    if (&rhs != this)
    {
      // Context.
      context_ = rhs.context_;
      
      // Deep copy of gates.
      gates_.clear();
      gates_.reserve(rhs.length());
      for (const auto& gate : rhs.gates_)  //  STYLE: Do we need const here?
      {
        gates_.push_back(gate->clone());
      }
      imprint_gates();  // Ensure gates get data from their new 'parent' circuit and not the old one.
    }

    marked_state = -5;  // Choosing an illegal value, rather than rhs.marked_state, catches bugs where a gate is...
                        // ...incorrectly 'imprinted'.
    return *this;
  }


  Circuit& Circuit::operator=(Circuit&& rhs)
  {
    // Check for self assignment
    if (&rhs != this)
    {
      context_ = std::exchange(rhs.context_, nullptr);  // Don't have to exchange with nullptr - the (default)
      gates_ = std::move(rhs.gates_);                   // ...destructor won't delete it anyway.
      imprint_gates();  // Ensure gates get data they need from their new 'parent' circuit and not the old one.
      marked_state = -6;  // Choosing an illegal value, rather than rhs.marked_state, catches bugs where a gate is
    }                     // ...incorrectly 'imprinted'.

    return *this;
  }


  bool Circuit::equivalentStructure(const Circuit& rhs) const
  {
    if (length() != rhs.length())
    {
      return false;
    }
    for (auto i = 0; i < length(); ++i)
    {
      if (!gates_[i]->equivalentStructure(*rhs.gates_[i]))
      {
        return false;
      }
    }
    return true;
  }


  bool Circuit::operator==(const Circuit& rhs) const
  {
    if (length() != rhs.length())
    {
      return false;
    }
    for (auto i = 0; i < length(); ++i)
    {
      if (*gates_[i] != *rhs.gates_[i])
      {
        return false;
      }
    }
    return true;
  }


  bool Circuit::operator!=(const Circuit& rhs) const
  {
    return !operator==(rhs);
  }


  bool Circuit::sortBefore(const Circuit& rhs) const
  {
    // To enable sorting, efficient counting of duplicates, etc. Must have a == b, sortBefore(a, b) or sortBefore(b, a).
    if (length() != rhs.length())
    {
      return length() < rhs.length();
    }
    for (auto i = 0; i < length(); ++i)
    {
      if (*gates_[i] != *rhs.gates_[i])
      {
        return gates_[i]->sortBefore(*rhs.gates_[i]);
      }
    }
    return false;
  }


  const int& Circuit::markedState() const
  {
    return marked_state;
  }


  void Circuit::setMarkedState(int markedState)
  {
    assert(0 <= markedState && markedState < dim(context_->numQbits()));
    marked_state = markedState;
  }


  int Circuit::length() const
  {
    return static_cast<int>(gates_.size());
  }


  int Circuit::numParameters() const
  {
    int count{0};
    for (auto& gate : gates_)
    {
      count += gate->numParameters();
    }

    return count;
  }


  const Gate& Circuit::gate(int pos) const
  {
    return *gates_[pos];
  }


  vector<double> Circuit::parameters() const
  {
    // Get all the angle parameters involved in the circuit.

    // Create the vector.
    vector<double> params;
    params.reserve(numParameters());

    // Populate it.
    for (auto& gate: gates_)
    {
      for (auto i = 0; i < gate->numParameters(); ++i)
      {
        params.push_back(gate->parameter(i));
      }
    }

    return params;
  }


  void Circuit::setParameters(const vector<double>& values)
  {
    // Set all the angle parameters involved in the circuit.
    if (values.size() != numParameters())
    {
      throw std::invalid_argument("Number of values provided does not match the number of parameters in Circuit::setParameters()");
    }
 
    auto i = 0;
    for (auto& gate : gates_)
    {
      for (auto j = 0; j < gate->numParameters(); ++j)
      {
        gate->setParameter(j, values[i++]);
      }
    }
  }


  void Circuit::random()
  {
    // Clear any existing gates.
    gates_.clear();

    // Determine how long the solution should be.
    auto length = randGeometric(context_->meanRandomLength());
    gates_.reserve(length);

    // Create a vector of random gates.
    for (auto i = 0; i < length; ++i)
    {
      // Get a random gate of random type and add it to the circuit.
      gates_.push_back(random_gate());  // random_gate() automatically 'imprints' the gate with its parent circuit.
    }
  }


  void Circuit::randomizeParameters()
  {
    for (auto& gate : gates_)
    {
      gate->randomizeParameters();
    }
  }


  bool Circuit::mutateReplaceGate()
  {
    // Simply switch out an existing gate for a randomly generated one.

    // First ensure that there is a gate to replace!
    if (length() == 0)
    {
      return false;  // At present, we don't mind if a call to one of these mutation functions results in no change.
    }

    // Select gate to remove.
    auto location = randInt(0, length());

    // Replace with a new randomly selected gate.
    gates_[location] = random_gate();  // random_gate() automatically 'imprints' the gate with its parent circuit.

    return true;
  }


  bool Circuit::mutateInsertGate()
  {
    // Insert a randomly generated gate at a random point in the circuit.

    // Select insertion location.
    auto location = randInt(0, length() + 1);  // A vector of 4 elements has 5 insertion points - before each...
                                               // ...element, and after the last.
    // Insert a random gate of random type.
    gates_.insert(gates_.begin() + location, random_gate());  // random_gate() automatically 'imprints' the gate with...
    return true;                                              // ...its parent circuit.
  }


  bool Circuit::mutateRemoveGate()
  {
    // Remove a randomly selected gate.

    // First ensure that there is a gate to remove!
    if (length() == 0)
    {
      return false;  // At present, we don't mind if a call to one of these mutation functions results in no change.
    }

    // Select gate to remove.
    auto location = randInt(0, length());

    // Remove it.
    gates_.erase(gates_.begin() + location);
    return true;
  }


  bool Circuit::mutateSwapGates()
  {
    // Select two existing gates in the circuit and switch their locations.

    // First ensure that there are at least two gates in the circuit.
    if (length() < 2)
    {
      return false;  // At present, we don't mind if a call to one of these mutation functions results in no change.
    }

    // Select the two gates
    auto location1 = randInt(0, length());
    int location2;
    do
    {
      location2 = randInt(0, length());
    }
    while (location1 == location2);

    // Swap them.
    std::swap(gates_[location1], gates_[location2]);
    return true;
  }


  bool Circuit::mutateMoveGate()
  {
    // Select an existing gate and move it to a new location.

    // First ensure that there are at least two gates in the circuit.
    if (length() < 2)
    {
      return false;  // At present, we don't mind if a call to one of these mutation functions results in no change.
    }

    // Select the gate to move and the position to move it to.
    auto mover = randInt(0, length());
    int newLocation;
    do
    {
      newLocation = randInt(0, length());
    }
    while (mover == newLocation);

    // Move it.
    auto b = gates_.begin();
    if (mover < newLocation)
    {
      std::rotate(b + mover, b + mover + 1, b + newLocation + 1);  // Last argument must point one beyond the...
    }                                                              // ...rotating section.
    else
    {
      std::rotate(b + newLocation, b + mover, b + mover + 1);
    }

    return true;
  }


  bool Circuit::mutateMutateGate()  // For 'adjust parameter and reoptimize' functionality.
  {
    // At present, this selects a gate at random, then adjusts a parameter (or more) if the gate has parameters. If, in
    // future, we wanted to ensure that this function actually changed the solution, we would have to modify it
    // somewhat. With the current overall method of performing mutations, this approach is fine for now.
    // (While other mutation operations, e.g. moveGate, swapGates, will make a change if they can (though even they
    // won't spot if they are swapping identical gates), this mutation operation simple gives up after one attempt.)

    // First ensure that there is a gate to mutate!
    if (length() == 0)
    {
      return false;  // At present, we don't mind if a call to one of these mutation functions results in no change.
    }

    // Select gate to mutate.
    auto location = randInt(0, length());
    auto& gate = gates_[location];  // Added & since 'gate' is now a unique_ptr.

    // Mutate
    if (gate->mutatable())
    {
      gate->mutate();
      return true;
    }

    return false;
  }


  State Circuit::simulate(const State& startState) const
  {
    // Simply apply each gate in turn to the provided startState, returning the result.
    State state(startState);
    for (auto& gate : gates_)
    {
      state = gate->applyTo(state);
    }

    return state;
  }


  vector<State> Circuit::simulateAndLog(const State& startState) const
  {
    // Like simulate(), but returns intermediate states as well. These are used to make later calculation of gradient
    // information more efficient.
    vector<State> visitedStates;
    visitedStates.reserve(length() + 1);

    visitedStates.push_back(startState);
    for (auto& gate : gates_)
    {
      visitedStates.push_back(gate->applyTo(visitedStates.back()));
    }

    return visitedStates;
  }


  vector<State> Circuit::simulateInverseAndLog(const State& startState) const
  {
    // Simulates the application of the inverse gates, in reverse order. Hence we can imagine the end state being
    // operated on by inverse gates until it turns into the start state. This is used to speed up calculation of the
    // gradient.
    //
    // Given target state t, start state s and gates ABCDEF, say we need to calculate t^FED'CBAs, where D' is the
    // derivative with respect to the D gate's angle parameter and t^ is the Hermitian conjugate. We can do this by
    // taking the stored result for CBAs, multiplying it by D' to get D'CBAs and then take the dot product with the
    // stored result for E^F^t. This means that, by doing both a forward and reverse simulation of the circuit and
    // storing the results, we can significantly reduce the effort required to calculate such derivatives and the
    // overall gradient.
    vector<State> visitedStates;
    visitedStates.reserve(length() + 1);

    visitedStates.push_back(startState);
    for (auto pos = static_cast<int>(length()) - 1; pos >= 0; --pos)
    {
      visitedStates.push_back(gates_[pos]->applyInvTo(visitedStates.back()));
    }

    return visitedStates;
  }


  cmplx Circuit::simulatedOverlap(const State& startState, const State& targetState) const
  {
    // Simulate the circuit starting with 'startState' and return the overlap of the end state with 'targetState'.
    // Function provided to assist problem implementers, since measures of circuit error will always be based on these
    // overlaps.
    auto endState = simulate(startState);
    return stateOverlap(endState, targetState);
  }


  std::pair<cmplx, vector<cmplx> > Circuit::simulatedOverlapAndGrad(const State& startState,
                                                                    const State& targetState) const
  {
    // Returns the overlap of the end state, resulting from simulating with the provided start state, and the
    // targetState. Also returns the gradient of this overlap.

    // Simulate the circuit both forwards and backwards, storing the intermediate results.
    auto visitedStatesForward = simulateAndLog(startState);
    auto visitedStatesBackward = simulateInverseAndLog(targetState);

    // To calculate the gradient, consider each gate parameter.
    vector<cmplx> gradOverlap(numParameters());
    int overallParamNum{0};
    for (int gateNum = 0; gateNum < length(); ++gateNum)
    {
      for (int gateParamNum = 0; gateParamNum < gate(gateNum).numParameters(); ++gateParamNum)
      {
        // Determine the overlap between, for example, E'DCBAs and F^G^t, which is t^GFE'DCBAs. Here s is the start
        // state, t is the target state and we are calculating the derivative with respect to some parameter of gate E
        // - hence the use of E'. F^ represents the inverse of F. Note how this uses the intermediate states from the
        // two simulations.
        gradOverlap[overallParamNum] =
             stateOverlap(gate(gateNum).applyGradTo(visitedStatesForward[gateNum], gateParamNum),
                          visitedStatesBackward[length() - gateNum - 1]);
        ++overallParamNum;
      }
    }

    // Determine the overlap itself.
    auto& endState = visitedStatesForward.back();
    cmplx overlap = stateOverlap(endState, targetState);

    return {overlap, gradOverlap};
  }


  void Circuit::makeCanonicalForm(int startPos)
  {
    // To avoid dealing with very many versions of the same circuit that merely have gates reordered, we convert to a
    // canonical form. Gates are reordered, as allowed given the rules encoded in swaps(), so that gates that are
    // 'sortedBefore' other gates come earlier where possible.
    //
    // Implemented as a simple modification to insertion sort.
    //
    // Parameter startPos allows us to save time - it indicates that gates before startPos are already in the correct
    // order, i.e. that insertion sort up to this position will simply leave the gates in place.

    // Create a vector of Gates, of the same length as the circuit, to act as storage for replacement gates, i.e. what
    // gates will become after some other gate is moved over them.
    GateSequence changedGates(length());

    // Consider each element in turn.
    for (int movingPos = startPos; movingPos < length(); ++movingPos)
    {
      auto movingGate = std::move(gates_[movingPos]);  // We std::move it back into the circuit's gate sequence later.
      auto bestPos = movingPos;             // Current best position to insert the gate and...
      auto bestGate = movingGate->clone();  // ...what the inserted gate will be.

      // Consider swapping with previous gates until this is no longer permitted.
      auto candidatePos = movingPos - 1;
      while (candidatePos >= 0)
      {
        if (gates_[candidatePos]->canSimplySwap(*movingGate))
        {
          changedGates[candidatePos] = {};  // Indicates that the gate does not change.
        }
        else if (gates_[candidatePos]->canSwapSwap(*movingGate))
        {
          changedGates[candidatePos] = gates_[candidatePos]->rightSwapMoverChange(*movingGate);
          auto changedMover = gates_[candidatePos]->leftSwapMoverChange(*movingGate);
          if (changedMover)
          {
            movingGate = std::move(changedMover);
          }
        }
        else if (gates_[candidatePos]->canHSwap(*movingGate))
        {
          changedGates[candidatePos] = gates_[candidatePos]->rightHMoverChange(*movingGate);
          auto changedMover = gates_[candidatePos]->leftHMoverChange(*movingGate);
          if (changedMover)
          {
            movingGate = std::move(changedMover);
          }
        }
        else if (gates_[candidatePos]->canRSwap(*movingGate))
        {
          changedGates[candidatePos] = gates_[candidatePos]->rightRMoverChange(*movingGate);
          auto changedMover = gates_[candidatePos]->leftRMoverChange(*movingGate);
          if (changedMover)
          {
            movingGate = std::move(changedMover);
          }
        }
        else
        {
          break;  // Cannot swap.
        }

        // We want to move the gate to the earliest position that currently contains a gate that is 'greater' than it,
        // under the constraint that it must be possible to swap the gate to that position.
        if (movingGate->sortBefore(*gates_[candidatePos]))
        {
          // Moving the gate to here results in a new circuit that 'sorts before' the current circuit (and previous
          // 'new' circuits).
          bestPos = candidatePos;
          bestGate = movingGate->clone();
        }

        --candidatePos;
      }

      // Having found the new position for the 'moving' gate, insert it and shift the jumped (and possibly modified)
      // gates up.
      int newMovingPos = movingPos;  // Needed if changes are made to jumped gates. Remember that 1 will be added to...
      for (auto changePos = movingPos; changePos > bestPos; --changePos)  // ...movingPos at loop end.
      {
        if (changedGates[changePos - 1])
        {
          gates_[changePos] = std::move(changedGates[changePos - 1]);

          // We may need to go back over bits that were previously sorted if any of these gates is changed.
          // (Note that if we were to make Hadamard::gate_type_id large, i.e. one less than the maximum (which is taken
          // by swap), then the following lines would be unnecessary - Hadamards would only jump back over Swaps, which
          // would be unchanged.)
          if (changePos - 1 == bestPos)
          {
            // The changed gate, being just after the insertion position, is correctly placed. The subsequent gate might
            // now wish to move back.
            newMovingPos = changePos;
          }
          else
          {
            // The changed gate might now wish to move back.
            newMovingPos = changePos - 1;
          }
        }
        else
        {
          gates_[changePos] = std::move(gates_[changePos - 1]);
        }
      }
      gates_[bestPos] = std::move(bestGate);
      movingPos = newMovingPos;
    }
  }


  void Circuit::makeUncanonical()
  {
    // While converting circuits to a canonical form has advantages for simplification and the logistics of cacheing, it
    // may result in a lack of diversity in the population for the genetic operators. Scrambling the solutions, using
    // legal swaps, should increase the number of possible circuits that could result from application of these
    // operators.

    vector<char> moved(length(), false);  // Array of bools indicating which gates have been moved.
                                          // (DO NOT USE vector<bool> here!!)
    // Consider each gate in turn
    int movingPos = 0;
    while (movingPos < length())
    {
      // Get range of permissible new places for the gate.
      int leftmost = leftmost_shift(movingPos, true);
      int rightmost = rightmost_shift(movingPos, true);

      // Select at random from this range
      int newPos = randInt(leftmost, rightmost + 1);

      // Perform the move, also rearranging the array indicating which positions have been considered.
      if (newPos > movingPos)
      {
        swap_gate_right(movingPos, newPos, true);
        std::rotate(moved.begin() + movingPos, moved.begin() + movingPos + 1, moved.begin() + newPos + 1);
      }
      else if (newPos < movingPos)
      {
        swap_gate_left(movingPos, newPos, true);
        std::rotate(moved.begin() + newPos, moved.begin() + movingPos, moved.begin() + movingPos + 1);
      }
      moved[newPos] = true;

      // Move to the next unmoved element.
      while (moved[movingPos] && movingPos < length())
      {
        ++movingPos;
      }
    }
  }


  shared_ptr<CachedStats> Circuit::simplifyCircuit(std::shared_ptr<CircuitCache> cache,
                                                   shared_ptr<CachedStats> emptyStats)
  {
    // My fifth circuit simplification routine.
    //
    // Perform steepest descent search until no 'simplification' can be found that either reduces cost or keeps cost
    // fixed while increasing the number of useful degrees of freedom for the numerical optimization. Having considered
    // all of the possible 'simplifications', the only ones that might realistically keep both cost and the number of
    // degrees of freedom fixed are swaps of gates such as XGates with Hadamards. We consider such changes to be merely
    // swaps when they do not change cost, rather than 'simplifications'. As a result, there is no need or inducement to
    // follow the steepest descent search with a full tree search of equivalent circuits.
    // NOTE: We assume that the circuit is in canonical form. In particular, all swaps should be at the end.
    //
    // At each step, we note the circuit visited (in canonical form) and add to the cache, with a pointer to a 'stats'
    // object that will eventually contain the simplified circuit and, after evaluation, the best fitness values
    // obtained, numerical optimization history, etc.
    //
    // This version of the simplification routine DOES NOT use the cache to shortcut simplification, as copying the
    // simplified circuit will usually result in an inequivalent circuit (with equivalent structure) due to the gate
    // parameters. (The exception, when the circuit has no parameterized gates, is handled by passing control to
    // simplifyStructure, which does allow shortcutting.)

    // If the circuit has no parameters, we can safely call simplifyStructure(), which might be quicker.
    if (numParameters() == 0)
    {
      return simplifyStructure(cache, emptyStats);
    }

    // Create a vector to store circuits visited during simplification, for adding later to the cache.
    vector<Circuit> visited;
    visited.reserve(100);  // MAGIC NUMBER: Probably excessive, though reallocation would be be expensive!
    visited.push_back(*this);  // Copying the circuit, but can't be avoided.

    auto bestReplacement = best_replacement();
    while (bestReplacement.valid())
    {
      // Make the simplification, making sure the circuit ends up in canonical form.
      apply_replacement(bestReplacement.leftPos(), bestReplacement.rightPos(), bestReplacement.meetPos());

      // Add adjusted circuit to the pile and look for further simplifications
      visited.push_back(*this);  // Copying the circuit, but can't be avoided.
      bestReplacement = best_replacement();
    }

    // No more simplifications can be made. Create a CachedStats object for the solution and add to the cache. (Note
    // that if the simplified circuit is already in the cache, cache.add() will notice this, adding only new links to
    // the existing CachedStats, rather than completing and linking to emptyStats.)
    if (cache)
    {
      return cache->add(visited, *this, emptyStats);
    }
    return emptyStats;  // Indicates that the solution is a new one (or wasn't in the cache.) (Must not return nullptr.)
  }


  shared_ptr<CachedStats> Circuit::simplifyStructure(std::shared_ptr<CircuitCache> cache,
                                                     shared_ptr<CachedStats> emptyStats)
  {
    // My fifth circuit simplification routine.
    //
    // Perform steepest descent search until no 'simplification' can be found that either reduces cost or keeps cost
    // fixed while increasing the number of useful degrees of freedom for the numerical optimization. Having considered
    // all of the possible 'simplifications', the only ones that might realistically keep both cost and the number of
    // degrees of freedom fixed are swaps of gates such as XGates with Hadamards. We consider such changes to be merely
    // swaps when they do not change cost, rather than 'simplifications. As a result, there is no need or inducement to
    // follow the steepest descent search with a full tree search of equivalent circuits.
    // NOTE: We assume that the circuit is in canonical form. In particular, all swaps should be at the end.
    //
    // At each step, we note the circuit visited (in canonical form) and add to the cache, with a pointer to a 'stats'
    // object that will eventually contain the simplified circuit and, after evaluation, the best fitness values
    // obtained, numerical optimization history, etc.
    //
    // This version of the simplification routine also uses the cache at each step to determine whether the circuit
    // structure has been seen before and, if it has, the simplified circuit is copied from the cache. Note that this
    // will mean that, while the simplified circuit structure is equivalant to the original, the circuit itself is
    // likely to differ due to different gate parameters. If we wish to create a simplified circuit that is equivalent
    // to the original, we call simplifyCircuit() which does not short-cut simplification in this way.

    // Have we already visited this circuit? If so, copy the simplified circuit from the cache.
    if (cache)
    {
      auto ptrToCachedStats = cache->find(*this);
      if (ptrToCachedStats)
      {
        *this = ptrToCachedStats->circuit();  // Do we need to copy other stuff too, e.g. objectives?
        return ptrToCachedStats;
      }
    }

    // This is a new circuit, or at least isn't in the cache. However, bear in mind that during simplification, we may
    // bump into a circuit that we have seen before and is in the cache.

    // Create a vector to store circuits visited during simplification, for adding later to the cache.
    vector<Circuit> visited;
    visited.reserve(100);  // MAGIC NUMBER: Probably excessive, though reallocation would be be expensive!
    visited.push_back(*this);  // Copying the circuit, but can't be avoided.

    auto bestReplacement = best_replacement();
    while (bestReplacement.valid())
    {
      // Make the simplification, making sure the circuit ends up in canonical form.
      apply_replacement(bestReplacement.leftPos(), bestReplacement.rightPos(), bestReplacement.meetPos());

      // Check to see if we have simplified to a circuit we have seen before.
      if (cache)
      {
        auto ptrToCachedStats = cache->find(*this);
        if (ptrToCachedStats)
        {
          *this = ptrToCachedStats->circuit();  // Do we need to copy other stuff too, e.g. objectives?
          cache->add(visited, ptrToCachedStats);  // Link visited circuits to the cached statistics for the...
          return ptrToCachedStats;                // ...simplified circuit.
        }
      }
      // This is still a circuit we haven't seen, or at least isn't stored in the cache. Add to the pile and look for
      // further simplifications
      visited.push_back(*this);  // Copying the circuit, but can't be avoided.
      bestReplacement = best_replacement();
    }

    // No more simplifications can be made. Since we have not seen the solution before, create a CachedStats object for
    // the solution and add to the cache.
    if (cache)
    {
      return cache->add(visited, *this, emptyStats);
    }
    return emptyStats;  // Indicates that the solution is a new one (or wasn't in the cache.) (Must not return nullptr.)
  }


  bool Circuit::reduce()
  {
    // Detect rotation gates that, due to the special value taken by the angle parameter, can be reduced to a cheaper
    // alternative or eliminated entirely. Perform this replacement/elimination. This is used after the application of
    // the numerical algorithm to optimize the gate parameters, since the numerical optimization may well determine that
    // a rotation gate is either not needed to produce the best results, or can be replaced by something like a ZGate.

    bool reduced = false;
    for (int pos = length() - 1; pos >= 0; --pos)  // Since reduction involves removing and inserting gates, it is...
    {                                              // ...easier to start from the end.
      if (gates_[pos]->canReduce())
      {
        GateSequence newSequence = gates_[pos]->reduction();
        auto replacementSize = static_cast<int>(newSequence.size());

        // Replace the old gate with the new ones.
        if (replacementSize == 0)
        {
          // If the gate is essentially the identity, we just remove it.
          gates_.erase(gates_.begin() + pos);
        }
        else
        {
          // Otherwise, we overwrite with the first/last (delete as appropriate!) gate in the new sequence and then
          // insert the rest.
          gates_[pos] = std::move(newSequence.back());
          gates_.insert(gates_.begin() + pos, make_move_iterator(newSequence.begin()),
                        make_move_iterator(newSequence.end() - 1));
        }

        reduced = true;
      }
    }

    // While we may be able to simplify, we leave this to the caller. However, we should still ensure that the circuit
    // ends up in canonical form.
    if (reduced)
    {
      makeCanonicalForm();
    }
    return reduced;
  }


  void Circuit::imprint_gates()
  {
    // Make sure each gate knows the identity of its 'parent' circuit. Used when circuits are copied (requiring the
    // copied gates to refer to the new parent) and created.
    for (auto& gate : gates_)
    {
      gate->imprint(*this);
    }
  }


  std::unique_ptr<Gate> Circuit::random_gate() const
  {
    auto newGate = context_->randomGate();
    newGate->imprint(*this);
    return newGate;
  }


  int Circuit::leftmost_shift(int pos, bool involveSwapGates) const
  {
    // Determine the leftmost position the gate can be pulled to using legal gate swaps. Parameter 'involveSwapGates'
    // indicates whether SwapSwaps are considered. Currently only used for making a circuit 'uncanonical'. Hence
    // efficiency is not a concern and 'involveSwapGates' is always 'true'.

    auto pulledPos = pos;
    auto movingGate = gates_[pos]->clone();
    while (pulledPos > 0)
    {
      if (gates_[pulledPos - 1]->canSimplySwap(*movingGate))
      {
        // Simple swap.
        --pulledPos;
      }
      else
      {
        // HSwaps, RSwaps and SwapSwaps may change the movingGate.
        unique_ptr<Gate> possibleChange;
        if (gates_[pulledPos - 1]->canHSwap(*movingGate))
        {
          // Can HSwap. Note possible change to the moving gate.
          possibleChange = gates_[pulledPos - 1]->leftHMoverChange(*movingGate);
        }
        else if (gates_[pulledPos - 1]->canRSwap(*movingGate))
        {
          // Can RSwap. Note possible change to the moving gate.
          possibleChange = gates_[pulledPos - 1]->leftRMoverChange(*movingGate);
        }
        else if (involveSwapGates && gates_[pulledPos - 1]->canSwapSwap(*movingGate))
        {
          // Can SwapSwap and we permit such operations. Note possible change to the moving gate.
          possibleChange = gates_[pulledPos - 1]->leftSwapMoverChange(*movingGate);
        }
        else
        {
          // No swap available.
          break;
        }

        // Can HSwap, RSwap or SwapSwap. Change movingGate if necessary and record gate identity in the movedGate array.
        if (possibleChange)  // If possibleChange is null, gate is unchanged.
        {
          movingGate = std::move(possibleChange);
        }
        --pulledPos;
      }
    }

    return pulledPos;
  }


  int Circuit::leftmost_shift(int pos, vector<const Gate*>& movedGate, vector<const Gate*>& pointers,
                              bool involveSwapGates) const
  {
    // Determine the leftmost position the gate can be pulled to using legal gate swaps. Also note how the moving gate
    // changes by recording the form of the gate in the movedGate array. Parameter 'involveSwapGates' indicates whether
    // SwapSwaps are considered. Only called, at present, when simplifying a circuit that is in canonical form. Hence
    // 'involveSwapGates' is always false. Efficiency may be a concern when the circuits being considered are
    // parameterless. We therefore keep gate creation to a minimum. These efficiency worries also explain the (risky)
    // use of bare pointers, rather than shared_ptrs, which were discovered to be too slow. The 'pointers' parameter
    // gets a copy of all the pointers that need to be deleted by the caller.

    auto pulledPos = pos;
    const Gate* movingGate = gates_[pos].get();
    movedGate[pos] = movingGate;
    while (pulledPos > 0)
    {
      if (gates_[pulledPos - 1]->canSimplySwap(*movingGate))
      {
        // Simple swap.
        --pulledPos;
        movedGate[pulledPos] = movingGate;
      }
      else
      {
        // HSwaps, RSwaps and SwapSwaps may change the movingGate.
        unique_ptr<Gate> possibleChange;
        if (gates_[pulledPos - 1]->canHSwap(*movingGate))
        {
          // Can HSwap. Note possible change to the moving gate.
          possibleChange = gates_[pulledPos - 1]->leftHMoverChange(*movingGate);
        }
        else if (gates_[pulledPos - 1]->canRSwap(*movingGate))
        {
          // Can RSwap. Note possible change to the moving gate.
          possibleChange = gates_[pulledPos - 1]->leftRMoverChange(*movingGate);
        }
        else if (involveSwapGates && gates_[pulledPos - 1]->canSwapSwap(*movingGate))
        {
          // Can SwapSwap and we permit such operations. Note possible change to the moving gate.
          possibleChange = gates_[pulledPos - 1]->leftSwapMoverChange(*movingGate);
        }
        else
        {
          // No swap available.
          break;
        }

        // Can HSwap, RSwap or SwapSwap. Change movingGate if necessary and record gate identity in the movedGate array.
        if (possibleChange)  // If possibleChange is null, gate is unchanged.
        {
          movingGate = possibleChange.release();  // Claims ownership. Dangerous - must delete somewhere, so we...
          pointers.push_back(movingGate);         // ...make a note of the pointer for future deletion.
        }
        --pulledPos;
        movedGate[pulledPos] = movingGate;
      }
    }

    return pulledPos;
  }


  int Circuit::rightmost_shift(int pos, bool involveSwapGates) const
  {
    // Determine the leftmost position the gate can be pulled to using legal gate swaps. Parameter 'involveSwapGates'
    // indicates whether SwapSwaps are considered. Currently only used for making a circuit 'uncanonical'. Hence
    // efficiency is not a concern and 'involveSwapGates' is always 'true'.

    auto pushedPos = pos;
    auto movingGate = gates_[pos]->clone();
    while (pushedPos < length() - 1)
    {
      if (movingGate->canSimplySwap(*gates_[pushedPos + 1]))
      {
        // Simple swap.
        ++pushedPos;
      }
      else
      {
        // HSwaps, RSwaps and SwapSwaps may change the movingGate.
        unique_ptr<Gate> possibleChange;
        if (movingGate->canHSwap(*gates_[pushedPos + 1]))
        {
          // Can HSwap. Note possible change to the moving gate.
          possibleChange = movingGate->rightHMoverChange(*gates_[pushedPos + 1]);
        }
        else if (movingGate->canRSwap(*gates_[pushedPos + 1]))
        {
          // Can RSwap. Note possible change to the moving gate.
          possibleChange = movingGate->rightRMoverChange(*gates_[pushedPos + 1]);
        }
        else if (involveSwapGates && movingGate->canSwapSwap(*gates_[pushedPos + 1]))
        {
          // Can SwapSwap and we permit such operations. Note possible change to the moving gate.
          possibleChange = movingGate->rightSwapMoverChange(*gates_[pushedPos + 1]);
        }
        else
        {
          // No swap available.
          break;
        }

        // Can HSwap, RSwap or SwapSwap. Change movingGate if necessary and record gate identity in the movedGate array.
        if (possibleChange)  // If possibleChange is null, gate is unchanged.
        {
          movingGate = std::move(possibleChange);
        }
        ++pushedPos;
      }
    }

    return pushedPos;
  }


  int Circuit::rightmost_shift(int pos, vector<const Gate*>& movedGate, vector<const Gate*>& pointers,
                               bool involveSwapGates) const
  {
    // Determine the rightmost position the gate can be pushed to using legal gate swaps. Also note how the moving gate
    // changes by recording the form of the gate in the movedGate array. Parameter 'involveSwapGates' indicates whether
    // SwapSwaps are considered. Only called, at present, when simplifying a circuit that is in canonical form. Hence
    // 'involveSwapGates' is always false. Efficiency may be a concern when the circuits being considered are
    // parameterless. We therefore keep gate creation to a minimum. These efficiency worries also explain the (risky)
    // use of bare pointers, rather than shared_ptrs, which were discovered to be too slow. The 'pointers' parameter
    // gets a copy of all the pointers that need to be deleted by the caller.

    auto pushedPos = pos;
    const Gate* movingGate = gates_[pos].get();
    movedGate[pos] = movingGate;
    while (pushedPos < length() - 1)
    {
      if (movingGate->canSimplySwap(*gates_[pushedPos + 1]))
      {
        // Simple swap.
        ++pushedPos;
        movedGate[pushedPos] = movingGate;
      }
      else
      {
        // HSwaps, RSwaps and SwapSwaps may change the movingGate.
        unique_ptr<Gate> possibleChange;
        if (movingGate->canHSwap(*gates_[pushedPos + 1]))
        {
          // Can HSwap. Note possible change to the moving gate.
          possibleChange = movingGate->rightHMoverChange(*gates_[pushedPos + 1]);
        }
        else if (movingGate->canRSwap(*gates_[pushedPos + 1]))
        {
          // Can RSwap. Note possible change to the moving gate.
          possibleChange = movingGate->rightRMoverChange(*gates_[pushedPos + 1]);
        }
        else if (involveSwapGates && movingGate->canSwapSwap(*gates_[pushedPos + 1]))
        {
          // Can SwapSwap and we permit such operations. Note possible change to the moving gate.
          possibleChange = movingGate->rightSwapMoverChange(*gates_[pushedPos + 1]);
        }
        else
        {
          // No swap available.
          break;
        }

        // Can HSwap, RSwap or SwapSwap. Change movingGate if necessary and record gate identity in the movedGate array.
        if (possibleChange)  // If possibleChange is null, gate is unchanged.
        {
          movingGate = possibleChange.release();  // Claims ownership. Dangerous - must delete somewhere, so we...
          pointers.push_back(movingGate);         // ...make a note of the pointer for future deletion
        }
        ++pushedPos;
        movedGate[pushedPos] = movingGate;
      }
    }

    return pushedPos;
  }


  Replacement Circuit::best_replacement() const
  {
    // Find the best replacement resulting from considering simplifications of pairs of gates, provided it improves cost
    // or keeps cost fixed but increases the number of useful degrees of freedom.
    // NOTE: We assume that the circuit is in canonical form. In particular, this means that all SwapGates should be at
    // the end (with the possible exception of Oracles).

    // First make a note of how far each gate may be shifted to the left and right and what the gate would look like in
    // each attainable position.
    vector<int> leftmost(length());
    vector<int> rightmost(length());
    vector<vector<const Gate*> > movedGate(length(), vector<const Gate*>(length(), nullptr));  // Try intrusive_ptrs?
    vector<const Gate*> pointers;
    pointers.reserve(2 * length());
    for (int pos = 0; pos < length(); ++pos)
    {
      leftmost[pos] = leftmost_shift(pos, movedGate[pos] , pointers, false);
      rightmost[pos] = rightmost_shift(pos, movedGate[pos] , pointers, false);
    }

    // Search through all possible replacements.
    long bestImprovement{0};
    bool bestExtraFreedom{false};
    Replacement bestReplacement;

    // Start by checking for any gates that can be moved to the start of the circuit and then discarded as redundant.
    for (int rightPos = 0; rightPos < length(); ++rightPos)
    {
      if (leftmost[rightPos] == 0)
      {
        // Gate can be shifted to the beginning.
        auto& newRightGate = movedGate[rightPos][0];
        if (newRightGate->cancelsAtStart())
        {
          auto improvement = newRightGate->cost();
          if (improvement > bestImprovement)
          {
            bestImprovement = improvement;
            bestReplacement = {-1, rightPos, -1};
          }
        }
      }
    }

    // Then look through all possible simplifications of gate pairs
    for (int leftPos = 0; leftPos < length() - 1; ++leftPos)
    {
      for (int rightPos = leftPos + 1; rightPos < length(); ++rightPos)
      {
        if (rightmost[leftPos] >= leftmost[rightPos] - 1)
        {
          // It is possible to move the two gates, using legal gate swaps, so that they meet, giving us the chance to
          // simplify.
          int meetPos = std::min(rightmost[leftPos], rightPos - 1);  // Choosing the rightmost position reduces work...
          auto& newLeftGate = movedGate[leftPos][meetPos];           // ...for makeCanonicalForm().
          auto& newRightGate = movedGate[rightPos][meetPos + 1];
          auto [improvement, extraFreedom] = newLeftGate->canSimplify(*newRightGate);
          if (improvement > bestImprovement ||
              improvement == bestImprovement && extraFreedom && !bestExtraFreedom)
          {
            bestImprovement = improvement;
            bestExtraFreedom = extraFreedom;
            bestReplacement = {leftPos, rightPos, meetPos};
          }
        }
      }
    }
    deletePointers(pointers);
    return bestReplacement;
  }


  void Circuit::swap_gate_right(int moverPos, int finalPos, bool involveSwapGates)
  {
    // Perform gate swaps such that the gate at 'moverPos' moves to the right until it reaches 'finalPos'. In the
    // process we adjust the moving gate and any gates jumped, if required. Parameter 'involveSwapGates' indicates
    // whether SwapSwaps are considered. (When simplifying a circuit that is in canonical form, they need not be
    // considered since all the SwapGates are out of the way. When 'un-canonicalizing' a circuit, SwapGates may be in
    // the way.)
    // This function assumes that it has already been determined that it is possible to move the gate into the desired
    // position. If not, it will throw an exception.

    assert(finalPos >= moverPos);
    auto mover = std::move(gates_[moverPos]);
    for (int pos = moverPos; pos < finalPos; ++pos)
    {
      if (mover->canSimplySwap(*gates_[pos + 1]))
      {
        // Just a simple swap
        gates_[pos] = std::move(gates_[pos + 1]);
      }
      else
      {
        // HSwaps, RSwaps and SwapSwaps may change the identify of one of the gates. (This code can handle both gates
        // changing, if necessary.)
        unique_ptr<Gate> newMover;
        unique_ptr<Gate> newJumpedGate;
        if (mover->canHSwap(*gates_[pos + 1]))
        {
          // We can HSwap. Make a note of gate changes.
          newMover = mover->rightHMoverChange(*gates_[pos + 1]);
          newJumpedGate = mover->leftHMoverChange(*gates_[pos + 1]);
        }
        else if (mover->canRSwap(*gates_[pos + 1]))
        {
          // We can RSwap. Make a note of gate changes.
          newMover = mover->rightRMoverChange(*gates_[pos + 1]);
          newJumpedGate = mover->leftRMoverChange(*gates_[pos + 1]);
        }
        else if (involveSwapGates && mover->canSwapSwap(*gates_[pos + 1]))
        {
          // We can SwapSwap, and permit such operations. Make a note of gate changes.
          newMover = mover->rightSwapMoverChange(*gates_[pos + 1]);
          newJumpedGate = mover->leftSwapMoverChange(*gates_[pos + 1]);
        }
        else
        {
          std::cerr << "Moving gate = " << *mover << std::endl;
          std::cerr << "Gate being jumped = " << *gates_[pos + 1] << std::endl;
          std::cerr << "Circuit:" << std::endl << *this;
          throw std::logic_error("Impossible to move gate to the desired position using the permitted swap operations");
        }

        // Perform the swap.
        if (newMover)
        {
          mover = std::move(newMover);
        }
        if (newJumpedGate)
        {
          gates_[pos] = std::move(newJumpedGate);
        }
        else
        {
          gates_[pos] = std::move(gates_[pos + 1]);
        }
      }
    }

    // Insert the moving gate into its final position.
    gates_[finalPos] = std::move(mover);
  }


  void Circuit::swap_gate_left(int moverPos, int finalPos, bool involveSwapGates)
  {
    // Perform gate swaps such that the gate at 'moverPos' moves to the left until it reaches 'finalPos'. In the process
    // we adjust the moving gate and any gates jumped, if required. Parameter 'involveSwapGates' indicates whether
    // SwapSwaps are considered. (When simplifying a circuit that is in canonical form, they need not be considered
    // since all the SwapGates are out of the way. When 'un-canonicalizing' a circuit, SwapGates may be in the way.)
    // This function assumes that it has already been determined that it is possible to move the gate into the desired
    // position. If not, it will throw an exception.

    assert(finalPos <= moverPos);
    auto mover = std::move(gates_[moverPos]);
    for (int pos = moverPos; pos > finalPos; --pos)
    {
      if (gates_[pos - 1]->canSimplySwap(*mover))
      {
        // Just a simple swap
        gates_[pos] = std::move(gates_[pos - 1]);
      }
      else
      {
        // HSwaps, RSwaps and SwapSwaps may change the identity of one of the gates. (This code can handle both gates
        // changing, if necessary.)
        unique_ptr<Gate> newMover;
        unique_ptr<Gate> newJumpedGate;
        if (gates_[pos - 1]->canHSwap(*mover))
        {
          // We can HSwap. Make a note of gate changes.
          newMover = gates_[pos - 1]->leftHMoverChange(*mover);
          newJumpedGate = gates_[pos - 1]->rightHMoverChange(*mover);
        }
        else if (gates_[pos - 1]->canRSwap(*mover))
        {
          // We can RSwap. Make a note of gate changes.
          newMover = gates_[pos - 1]->leftRMoverChange(*mover);
          newJumpedGate = gates_[pos - 1]->rightRMoverChange(*mover);
        }
        else if (involveSwapGates && gates_[pos - 1]->canSwapSwap(*mover))
        {
          // We can SwapSwap, and permit such operations. Make a note of the gate changes.
          newMover = gates_[pos - 1]->leftSwapMoverChange(*mover);
          newJumpedGate = gates_[pos - 1]->rightSwapMoverChange(*mover);
        }
        else
        {
          throw std::logic_error("Impossible to move gate to the desired position using the permitted swap operations");
        }

        // Perform the swap.
        if (newMover)
        {
          mover = std::move(newMover);
        }
        if (newJumpedGate)
        {
          gates_[pos] = std::move(newJumpedGate);
        }
        else
        {
          gates_[pos] = std::move(gates_[pos - 1]);
        }
      }
    }

    // Insert the moving gate into its final position.
    gates_[finalPos] = std::move(mover);
  }


  void Circuit::apply_replacement(int leftPos, int rightPos, int meetPos)
  {
    // The gates at leftPos and rightPos can be moved so as to meet at meetPos and a simplication/replacement can be
    // performed. Gates will change in the process of moving if HSwaps or RSwaps are used.
    // If leftPos (and meetPos) are -1, this means that the gate at rightPos is moved to the start and then eliminated
    // instead.

    if (leftPos == -1)
    {
      // Move the right gate to the start, adjusting any gates jumped if a HSwap or RSwap is used, then remove the gate.
      swap_gate_left(rightPos, 0, false);
      gates_.erase(gates_.begin());
    }
    else
    {
      // Move the left gate to the right, adjusting it and any gates jumped if an HSwap or RSwap is used.
      swap_gate_right(leftPos, meetPos, false);

      // Move the right gate to the left, adjusting it and any gates jumped if a HSwap or RSwap is used.
      swap_gate_left(rightPos, meetPos + 1, false);

      // Now get the replacement sequence.
      GateSequence newSequence = gates_[meetPos]->simplification(*gates_[meetPos + 1]);
      auto replacementSize = static_cast<int>(newSequence.size());
      assert(replacementSize < 3);  // Circuit shouldn't grow. This means we can avoid using expensive insert() on...
                                    // ...the circuit's gate sequence.

      // Replace the old gates with the new ones.
      gates_.erase(gates_.begin() + meetPos + replacementSize, gates_.begin() + meetPos + 2);
      std::copy(make_move_iterator(newSequence.begin()), make_move_iterator(newSequence.end()),
                gates_.begin() + meetPos);
    }

    // Put the new circuit into canonical form.
    // (Should we do this here, or outside? Since the circuit is in canonical form at the start of this function, I
    // decided that this tidying up belonged here, at least for now.)
    makeCanonicalForm();
  }


  size_t Circuit::hash() const
  {
    size_t h{0};
    for (const auto& gate : gates_)
    {
      h *= context_->numPermittedGateOptions();
      h += gate->gateOptionId();
    }

    return h;
  }


  bool Circuit::gates_available() const
  {
    // A simple check on the availability of all the gates in the circuit.
    for (auto& gate : gates_)
    {
      if (!gate->available())
      {
        cerr << "Error: " << *gate << "is unavailable, yet appears in the circuit after simplification." << endl;
        return false;
      }
    }
    return true;
  }


  void Circuit::output(ostream& out) const
  {
    for (auto& gate : gates_)  // Must have & as 'gate' is now a unique_ptr.
    {
      if (gate)
      {
        out << *gate;
      }
      else
      {
        out << "Null gate!" << endl;
      }
    }
  }


  void crossover(Circuit& lhs, Circuit& rhs)
  {
    // At present, this just performs a two point crossover.

    int lhsCross1 = randInt(0, lhs.length() + 1);
    int lhsCross2 = randInt(0, lhs.length() + 1);
    if (lhsCross1 > lhsCross2)
    {
      std::swap(lhsCross1, lhsCross2);
    }

    int rhsCross1 = randInt(0, rhs.length() + 1);
    int rhsCross2 = randInt(0, rhs.length() + 1);
    if (rhsCross1 > rhsCross2)
    {
      std::swap(rhsCross1, rhsCross2);
    }

    GateSequence newLhs;
    newLhs.insert(newLhs.end(), make_move_iterator(lhs.gates_.begin()),
                  make_move_iterator(lhs.gates_.begin() + lhsCross1));
    newLhs.insert(newLhs.end(), make_move_iterator(rhs.gates_.begin() + rhsCross1),
                  make_move_iterator(rhs.gates_.begin() + rhsCross2));
    newLhs.insert(newLhs.end(), make_move_iterator(lhs.gates_.begin() + lhsCross2),
                  make_move_iterator(lhs.gates_.end()));

    GateSequence newRhs;
    newRhs.insert(newRhs.end(), make_move_iterator(rhs.gates_.begin()),
                  make_move_iterator(rhs.gates_.begin() + rhsCross1));
    newRhs.insert(newRhs.end(), make_move_iterator(lhs.gates_.begin() + lhsCross1),
                  make_move_iterator(lhs.gates_.begin() + lhsCross2));
    newRhs.insert(newRhs.end(), make_move_iterator(rhs.gates_.begin() + rhsCross2),
                  make_move_iterator(rhs.gates_.end()));

    lhs.gates_ = std::move(newLhs);
    rhs.gates_ = std::move(newRhs);

    // Ensure that any gates that have moved are directed to their new 'parent' for any circuit data (e.g. marked state)
    // that they might need.
    lhs.imprint_gates();
    rhs.imprint_gates();
  }


  ostream& operator<<(ostream& out, const Circuit& circuit)
  {
    circuit.output(out);
    return out;
  }

  //----------------------------------------------------------------------------------------------------------------

  // A base class for all quantum circuit Solution classes. Doesn't need to implement everything required by the genetic
  // algorithms, as classes such as NSGA2Solution will be templated on the derived Solution classes, not this base
  // class. (Indeed, the genetic algorithm code has no need for a Solution hierarchy, being template based.)

  Solution::Solution(const Problem& problem) :
  problem_(problem),
  circuit_(problem.circuitContext()),  // Member 'evaluated_' set to false by in-class initializer.
  error_(problem.numObj() - 1, -1.0),
  gradient_(numParameters(), 0.0),  // It feels risky assigning a reasonable looking gradient here, when the gradient...
  cache_(problem.cache()),          // ...should be 'unassigned'.
  stats_(nullptr),
  front_(problem.front())
  {
  }


  bool Solution::operator==(const Solution& rhs) const
  {
    // Equality operator just considers the solution itself, i.e. the sequence of gates. This means that we can use this
    // to detect that an unevaluated solution is, in fact, the same as an evaluated one.

    return circuit_ == rhs.circuit_;
  }


  bool Solution::operator!=(const Solution& rhs) const
  {
    return !operator==(rhs);
  }


  bool Solution::sortBefore(const Solution& rhs) const
  {
    // To enable sorting, efficient counting of duplicates, etc. Must have a == b, sortBefore(a, b) or sortBefore(b, a).

    return circuit_.sortBefore(rhs.circuit_);
  }


  void Solution::random()
  {
    // Create a random circuit and mark as unevaluated and put it into canonical form.
    circuit_.random();
    evaluated_ = false;
    circuit_.makeCanonicalForm();
  }


  void Solution::mutate(double mutateProb)
  {
    // The parameter 'mutateProb' was initially designed (in the MOMH code) to mean the probability of mutating each
    // individual bit of the solution. That isn't appropriate here, as there are many different forms of mutation: gate
    // replacement, gate mutation, gate insertion, gate removal, etc. For now, we will just use mutateProb to mean the
    // chance of applying each of the operators.

    // The solution is marked as unevaluated by each of the individual mutation operations, if they result in a change
    // to the solution. Also, the circuit involved is first scrambled using swaps that do not change its overall action,
    // using makeUncanonical, and unscrambled afterwards using makeCanonicalForm(). This increases the range of circuits
    // that can be created via mutation and crossover.

    circuit_.makeUncanonical();  // Scramble (without changing circuit output) to increase diversity.
    if (rand01() < mutateProb)
    {
      if (circuit_.mutateReplaceGate())
        evaluated_ = false;
    }
    if (rand01() < mutateProb)
    {
      if (circuit_.mutateInsertGate())
        evaluated_ = false;
    }
    if (rand01() < mutateProb)
    {
      if (circuit_.mutateRemoveGate())
        evaluated_ = false;
    }
    if (rand01() < mutateProb)
    {
      if (circuit_.mutateSwapGates())
        evaluated_ = false;
    }
    if (rand01() < mutateProb)
    {
      if (circuit_.mutateMoveGate())
        evaluated_ = false;
    }
    if (rand01() < mutateProb)  // Added to enable 'change parameter and reoptimize' functionality.
    {
      if (circuit_.mutateMutateGate())
        evaluated_ = false;
    }
    circuit_.makeCanonicalForm();
  }


  int Solution::numParameters() const
  {
    return circuit_.numParameters();
  }


  vector<double> Solution::parameters() const
  {
    return circuit_.parameters();
  }


  void Solution::setParameters(const vector<double> &values)
  {
    // Note that the need to set 'evaluated_' to false makes providing this function a better option than simply
    // providing access to 'circuit_'.
    circuit_.setParameters(values);
    evaluated_ = false;
  }


  shared_ptr<CachedStats> mostReduced(shared_ptr<CachedStats> stats)
  {
    // Thread-safe function to get the most reduced circuit, by following links in the cache. I assume that the locking
    // is necessary, i.e. that it would be possible to attempt to read reduced() while it is being written to, resulting
    // in a scrambled mess, otherwise. Note that we require that any 'writer's lock' that this thread has on stats be
    // released.
    auto rLock = stats->getReadersLock();
    while (stats->reduced())
    {
      stats = stats->reduced();
      rLock = stats->getReadersLock();
    }

    return stats;
  }


  void Solution::evaluateObjectives()
  {
    // Function required by the MOMH library. This evaluates the objectives for the circuit structure, not just the
    // circuit. Hence it includes simplification, numerical optimization, reduction etc.

    // First check whether this circuit is already evaluated and optimized - no point repeating the task!
    if (evaluated_)  // Warning: We assume that evaluated => simplified and optimized.
    {
      return;  // Don't waste time evaluating already evaluated solutions!
    }

    // Simplify the circuit. Circuits visited during simplification are added to the cache. The cache is also used to
    // spot if this circuit structure, or any resulting from the simplification, has previously been simplified and
    // evaluated. If so, we can simply look up the simplified circuit and error values.
    auto newStats = make_shared<CachedStats>();  // Used if the simplified circuit is new.
    auto newWLock = newStats->getWritersLock();  // Only this thread gets access until it is filled.
    stats_ = circuit_.simplifyStructure(cache_, newStats);
    assert(stats_);

    evaluate_cost();  // Counting up gate cost is quick.

    if (stats_ != newStats)
    {
      // The simplified circuit has already been visited. It has probably also been evaluated. If not, some other thread
      // is evaluating it (or about to), and it makes sense for this thread to wait until this is completed. While, in
      // future, we may have options to perform additional numerical optimization runs on the circuit, for now we assume
      // that the previous numerical optimization performed perfectly and check whether the circuit was also 'reduced'.
      // If so, we copy across the reduced circuit and the error values.

      // Check for cached reductions of the previously visited circuit.
      stats_ = mostReduced(stats_);
      auto statsRLock = stats_->getReadersLock();
      circuit_ = stats_->circuit();

      // Copy across the circuit error data.
      evaluated_ = true;
      error_ = stats_->error();

      return;
    }

    // Simplified circuit is new.

    // Rename (or rather transfer) lock from newWLock to statsWLock. (stats_ and newStats are the same.)
    auto statsWLock = std::move(newWLock);

    // Don't evaluate very big solutions. Also don't evaluate circuits if a good enough solution has already been found
    // at that circuit's cost.
    if (circuit_.length() > problem_.maxLength() ||  // Just too big!!
        front_.best(cost() - problem_.expectedReduction()) < problem_.targetError())  // Good enough found already at...
    {                                                                                 // ...lesser cost.
      // If the solution isn't simply 'too big', then we must have already found a solution of this size which is good
      // enough. In either case we may wish to focus the search on smaller circuits. For this circuit, this means giving
      // it a bad error without evaluating or performing numerical optimization.
      // This runs the risk of missing out on one chance of producing a good small circuit - the possibility of reducing
      // a larger circuit - if the estimate of the expected reduction is too small.
      // (Possible alternative: consider applying mutateRemoveGate() until the circuit is small enough. (Use
      // simplification at each step?)
      assign_worst_error();
      stats_->setEvaluated(error_);
      stats_->setReduced(nullptr);
      evaluated_ = true;
      return;
    }

    // If there are gate parameters, perform numerical optimization and then check for gate 'reductions'.
    if (problem_.randomizeParameters())
    {
      // Given a new circuit, randomize the gate parameters. Since these are likely to be near optimal for the parent
      // circuit, they may be expected to be good for this circuit too. Hence randomizing them tends to result in slower
      // optimization, but better exploration of parameter space.
      circuit_.randomizeParameters();
    }
    optimize_and_cache();  // Numerically optimize (if there are gate parameters), evaluate and cache, both in...
                           // stats_ and front_.
    // Having performed any required numerical optimization, look to see whether the angles of any of the parameterized
    // gates has become (approximately) zero or some other special value that allows replacement by a cheaper
    // alternative.
    while (circuit_.reduce())
    {
      // Circuit has been 'reduced', so now it may also be possible to further simplify the circuit. Circuits occuring
      // during simplification are added to the cache, with links to a NEW CachedStats object. (We want to keep the
      // numerical optimization data stored in the original CachedStats, in case we choose to re-optimize.) The old
      // CachedStats object also gets a link to the new CachedStats for the reduced circuit. We then redirect 'stats_'
      // to point to the new CachedStats object.
      auto newReducedStats = make_shared<CachedStats>();  // Used if the simplified reduced circuit is new.
      auto newReducedWLock = newReducedStats->getWritersLock();  // Only this thread gets access to the new...
                                                                 // ...CachedStats until it is filled.
      auto reducedStats = circuit_.simplifyCircuit(cache_, newReducedStats);  // Simplify reduced circuit and get...
                                                                              // ...pointer to the CachedStats.
      assert(reducedStats);
      stats_->setReduced(reducedStats);  // Link old CachedStats to the new. This is the last edit to the CachedStats...
      statsWLock.unlock();               // ...pointed to by stats_, so we can unlock it and redirect stats to point...
      stats_ = reducedStats;             // ...to the CachedStats for the reduced solution.

      evaluate_cost();

      if (stats_ != newReducedStats)
      {
        // The reduced circuit has already been visited.
        statsWLock = stats_->getWritersLock();
        if (primary_error() < stats_->primaryError())
        {
          // Circuit after numerical optimization (but before reduction) has better overall error than that stored in
          // the cache. Briefly run numerical optimization to fix any minor loss of quality during reduction) and then
          // update cache and front. (Gate parameters should already be close to a local minimum.)
          optimize_and_cache();
        }
        else
        {
          // Overall error of cached solution is better. Grab error values from the cache. Note, however, that the
          // cached data will involve different gate parameter values which may have lead to further reductions. So we
          // ensure that we copy across the most reduced version of the stored circuit first.
          statsWLock.unlock();  // We will not be writing to the CachedStats object, and do not mind if it gets...
          stats_ = mostReduced(stats_);  // ...updated by another thread.
          auto statsRLock = stats_->getReadersLock();
          circuit_ = stats_->circuit();

          // Copy across the circuit error data
          evaluated_ = true;
          error_ = stats_->error();
          break;  // Have copied across the better, reduced circuit from the cache: there can be no further reductions.
        }
      }
      else
      {
        // Circuit after reduction is new. Briefly run numerical optimization, to fix any minor loss of quality during
        // reduction, and then update cache and front. (Gate parameters should already be close to a local minimum.)
        statsWLock = std::move(newReducedWLock);
        optimize_and_cache();
      }
    }
  }


  double Solution::objective(int objNum) const
  {
    assert(evaluated_);
    assert(0 <= objNum && objNum < problem_.numObj());
    if (objNum < problem_.numObj() - 1)
    {
      return error_[objNum];
    }
    return total_cost;  // Casts total_cost to double.
  }


  long Solution::cost() const
  {
    return total_cost;
  }


  void Solution::output(ostream& out) const
  {
    out << "Circuit:" << endl << circuit_;

    if (evaluated_)
    {
      out << "Error: " << error_[0];
      for (auto i = 1; i < problem_.numObj() - 1; ++i)
      {
        out << ", " << error_[i];
      }
      out << endl;
      out << "Gate cost: " << total_cost << endl;
    }
    else
    {
      out << "Solution unevaluated." << endl;
    }
  }


  void Solution::outputGradient(std::ostream &out) const
  {
    // STYLE: Use an algorithm?
    if (numParameters() == 0)
    {
      out << "Circuit has no gate parameters to optimize." << std::endl;
    }
    else
    {
      out << "Gradient: " << gradient_[0];
      for (auto i = 1; i < numParameters(); ++i)
      {
        out << ", " << gradient_[i];
      }
      out << std::endl;
    }
  }


  void Solution::evaluate_cost()
  {
    // The easy bit of evaluating a circuit - adding up the cost of the gates.

    // We could use 'accumulate', as follows, but quite frankly I think the manual code is clearer.
    // total_cost = std::accumulate(gates_.begin(), gates_.end(), 0,
    //                           [*this](int sum, auto& element){auto gate = element.first; return sum + gate.cost();});
    total_cost = 0;
    for (auto i = 0; i < circuit_.length(); ++i)
    {
      total_cost += circuit_.gate(i).cost();
    }
  }


  double Solution::primary_error() const
  {
    // The primary error function, i.e. that minimized by the numerical optimizer, is the first in the list.
    return error_[0];
  }


  void Solution::optimize_parameters()
  {
    // Optimize the gate parameters using the LBFGS++ implementation of L-BFGS.

    // Create solver object.
    LBFGSParam<double> param;
    param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
    param.m = std::min(numParameters(), 40);
    param.epsilon = 0.000001;
    //    param.max_linesearch = 200;
    param.targetTest = true;
    param.past = 3;
    param.target = front_.best(cost());
    param.goodEnoughTest = true;
    param.goodEnough = problem_.targetError();
    param.iterationQuota = problem_.iterationQuota();
    LBFGSSolver<double> solver(param);

    // Set up objective (and gradient) function object and initial solution.
    OverallErrorAndGradient objFun(*this);
    vector<double> stdParameters = parameters();                   // vectorToVectorXd needs a non-const...
    VectorXd currentParameters = vectorToVectorXd(stdParameters);  // ...(hence non-temporary) std::vector to link to.
    double fx;

    // Run LBFGS
    solver.minimize(objFun, currentParameters, fx);

    // Check for hiccups.
    for (int i = 0; i < numParameters(); ++i)
    {
      if (currentParameters[i] != currentParameters[i])
      {
        cerr << "One of the parameters is no longer a number!!!" << endl;
        exit(EXIT_FAILURE);
      }
    }
    setParameters(vectorXdToVector(currentParameters));
  }


  void Solution::optimize_and_cache()
  {
    // Perform the numerical optimimization and the cacheing of the results. If there are no gate parameters, simply
    // evaluate and cache.

    // First, numerical optimization, if there are gate paremeters, cacheing the best parameters found.
    if (numParameters() > 0)
    {
      optimize_parameters();
      if (cache_)
      {
        stats_->setParameters(parameters());
      }
    }

    // Re-evaluate. (The numerical optimization (if run) has, of course already calculated the overall and worst case
    // errors. However the last set of parameters evaluated may not have been the best. While the LBFGS object can give
    // us the overall error, it doesn't have access to the worst case.)
    evaluate_error();

    // Update cached statistics and the front.
    if (cache_)
    {
      stats_->setEvaluated(error_);
      stats_->setReduced(nullptr);
    }
    front_.add(cost(), primary_error());
    evaluated_ = true;
  }


  void crossover(Solution& lhs, Solution& rhs)
  {
    // Scramble circuits using swaps that do not change circuit output. This may help increase diversity.
    lhs.circuit_.makeUncanonical();
    rhs.circuit_.makeUncanonical();

    // Perform the crossover.
    crossover(lhs.circuit_, rhs.circuit_);
    lhs.evaluated_ = false;
    rhs.evaluated_ = false;

    // Unscramble, i.e. put circuits in canonical form.
    lhs.circuit_.makeCanonicalForm();
    rhs.circuit_.makeCanonicalForm();
  }


  ostream& operator<<(ostream& out, const Solution& solution)
  {
    solution.output(out);
    return out;
  }


}
