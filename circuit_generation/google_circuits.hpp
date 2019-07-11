#include "ps_stabilizer/core.hpp"
#include "ps_stabilizer/chstabilizer.hpp"
#include "ps_stabilizer/rng.hpp"

#include "nlohmann/json.hpp"

#define _USE_MATH_DEFINES

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

enum class GateEnum {I, H, CZ, SX, SY, Rot};

std::unordered_map<GateEnum,std::string> gate_strs = {
  {GateEnum::H, "H"},
  {GateEnum::SX, "SX"},
  {GateEnum::SY, "SY"}
};

struct Gate;

using uint_t = StabilizerSimulator::uint_t;
using pauli_t = StabilizerSimulator::pauli_t;
using complex_t = StabilizerSimulator::complex_t;
using json_t = nlohmann::json;
using branch_t = std::pair<double, pauli_t>;
using layer_t = std::vector<GateEnum>;
using cz_t = std::pair<unsigned, unsigned>;
using gate_ptr = std::shared_ptr<Gate>;

layer_t EMPTY = {};

pauli_t IDENTITY;

static uint_t zer = 0ULL;
static uint_t one = 1ULL;

static const complex_t pi_over_8_phase(0., M_PI/8);
static const complex_t root_omega = std::exp(pi_over_8_phase);
static const complex_t root_omega_star = std::conj(root_omega);

struct Gate
{
    unsigned qubit;
    unsigned target;
    GateEnum type;

    Gate() = default;
    Gate(unsigned qubit_, GateEnum type_) : qubit(qubit_), target(0), type(type_) {};
    Gate(unsigned qubit_, unsigned target_): qubit(qubit_), target(target_), type(GateEnum::CZ) {}; 
    Gate(const Gate& other) : qubit(other.qubit), target(other.target), type(other.type) {};
    virtual ~Gate() = default;
    virtual void commute_past(gate_ptr other) { return; };

    operator std::string() const
    {
      switch(type)
      {
        case GateEnum::CZ:
        {
          std::string ctrl = std::to_string(qubit);
          std::string trgt = std::to_string(target);
          return "CZ["+ctrl+","+trgt+"];";
        }
        case GateEnum::I:
        {
          return "";
        }
        case GateEnum::Rot:
        {
          throw std::runtime_error("Case handled by the derived class.");
          break;
        }
        default:
        {
          std::string gstring = gate_strs[type];
          return gstring + "["+std::to_string(qubit)+"];";
        }
      }
    };

    virtual std::string as_qasm()
    {
      switch (type)
      {
        case GateEnum::H:
          return "h q0[" + std::to_string(qubit) + "];\n";
          break;
        case GateEnum::SX:
        {
          std::string qstring = std::to_string(qubit);
          return "h q0[" + qstring + "];\ns q0[" + qstring + "];\nh q0[" + qstring + "];\n";
          break;
        }
        case GateEnum::SY:
        {
          std::string qstring = std::to_string(qubit);
          return "s q0[" + qstring + "];\nh q0[" + qstring + "];\ns q0[" + qstring + "];\nh q0[" + qstring + "];\nsdg q0[" + qstring + "];\n";
          break;
        }
        case GateEnum::CZ:
        {
          return "cz q0[" + std::to_string(qubit) + ",q0[" + std::to_string(target) +"];\n";
          break;
        }
        default:
          return "";
          break;
      }
    }
};

struct Rotation : public Gate 
{
  pauli_t rotation_operator;
  double rotation_angle;
  double branch_angle; 

  double p_threshold;
  std::array<complex_t, 2> coeffs;

  Rotation() : rotation_angle(M_PI/8.),
  branch_angle(M_PI/4.), p_threshold(0.5)
  {
    type = GateEnum::Rot;
  };
  
  Rotation(unsigned qubit_) : Rotation() 
  {
    qubit = qubit_;
    rotation_operator.X = zer;
    rotation_operator.e = 3;
    rotation_operator.Z = (one << qubit_);
  };

  Rotation(const Rotation& other) : Gate(other)
  {
    rotation_operator = other.rotation_operator;
    rotation_angle = other.rotation_angle;
    branch_angle = other.branch_angle;
    p_threshold = other.p_threshold;
    coeffs = other.coeffs;
  }

  ~Rotation() = default;

  void commute_past(gate_ptr other) override
  {
    switch(other->type)
    {
      case GateEnum::H:
      {
        uint_t shift = (one << other->qubit);
        uint_t x_temp = (rotation_operator.X & (~shift)) ^ (rotation_operator.Z & shift);
        uint_t z_temp = (rotation_operator.Z & (~shift)) ^ (rotation_operator.X & shift);
        if (!!(rotation_operator.X & rotation_operator.Z & shift))
        {
          rotation_operator.e ^= 2U;
        }
        rotation_operator.X = x_temp;
        rotation_operator.Z = z_temp;
        break;
      }
      case GateEnum::SX:
      {
        uint_t shift = (one << other->qubit);
        if (!!(rotation_operator.Z & shift)) //Only update Z/Y
        {
          rotation_operator.e = (rotation_operator.e + 3)%4;
          rotation_operator.X ^= shift;
        }
        break;
      }
      case GateEnum::SY:
      {
        uint_t shift = (one << other->qubit);
        if(!!((rotation_operator.X ^ rotation_operator.Z) & shift)) //Only update if X or Z
        {
          rotation_operator.e ^= (!!(rotation_operator.X & shift)) * 2U;
          rotation_operator.X ^= shift;
          rotation_operator.Z ^= shift;
        }
        break;
      }
      case GateEnum::CZ:
      {
        uint_t control_shift = (one << other->qubit);
        uint_t target_shift = (one << other->target);
        //Update only if we act as X or Y on one of the qubits.
        if(!!(rotation_operator.X & control_shift) || !!(rotation_operator.X & target_shift))
        {
          rotation_operator.Z ^= (!!(rotation_operator.X & control_shift)) * target_shift;
          rotation_operator.Z ^= (!!(rotation_operator.X & target_shift)) * control_shift;
          rotation_operator.e ^= (!!(rotation_operator.X & control_shift) && !!(rotation_operator.X & target_shift))*2U;
        }
        break;
      }
      case GateEnum::I:
      {
        break;
      }
      default:
      {
        throw std::runtime_error("Commutation should stop rather than reach this point...");
        break;
      }
    }
  };

  std::string to_string(unsigned n_qubits)
  {
    if(rotation_operator.weight() == 1)
    {
      std::string out = "";
      unsigned pauli = !!(rotation_operator.X & (one << qubit))*2 + !!(rotation_operator.Z & (one << qubit))*1;
      switch (pauli)
      {
        case 0:
          out += "I";
          break;
        case 1:
          out += "Z";
          break;
        case 2:
          out += "X";
          break;
        case 3:
          out += "Y";
          break;
      }
      out += "["+std::to_string(qubit) + "@" + std::to_string(rotation_angle) + "];";
      return out;
    }
    std::string out = rotation_operator.to_string(n_qubits);
    out += "@" + std::to_string(rotation_angle) + ";";
    return out;
  };

  std::string as_qasm() override //Only implemented for uncompiled circuits
  {
    return "u1(" + std::to_string(rotation_angle) + ") q0[" + std::to_string(qubit) + "];\n";
  }
};

void clifford_recompilation(std::vector<gate_ptr>& gate_sequence)
{
  std::vector<unsigned> non_clifford_positions;
  unsigned n_gates = gate_sequence.size();
  for(unsigned i = n_gates-1; i<0; i--)
  {
    if(gate_sequence[i]->type == GateEnum::Rot)
    {
      non_clifford_positions.push_back(i);
    }
  }
  for(auto position : non_clifford_positions)
  {
    auto non_clifford = gate_sequence[position];
    unsigned current_pos = position;
    unsigned next_pos = position+1;
    auto next = gate_sequence[next_pos];
    while(next_pos < gate_sequence.size())
    {
      if(next->type == GateEnum::H || next->type == GateEnum::Rot)
      {
        break;
      }
      non_clifford->commute_past(next);
      gate_sequence[current_pos] = next;
      current_pos++;
      next_pos++;
    }
    gate_sequence[current_pos] = non_clifford;
  }
}

layer_t get_cz_layer(unsigned n_qubits, unsigned layer_index, std::vector<std::vector<cz_t>>& cz_schema, std::vector<gate_ptr>& gates)
{
  layer_t layer(n_qubits, GateEnum::I);
  for(auto cz : cz_schema[layer_index])
  {
    layer[cz.first] = GateEnum::CZ;
    layer[cz.second] = GateEnum::CZ;
    gates.push_back(std::make_shared<Gate>(cz.first, cz.second));
  }
  return layer;
}

layer_t get_cz_layer(unsigned n_qubits, double cz_fraction, std::vector<gate_ptr>& gates, layer_t& last_layer=EMPTY)
{
  layer_t layer(n_qubits, GateEnum::I);
  unsigned cz_pairs = std::rint(0.5*cz_fraction*n_qubits);
  unsigned pairs = 0;
  while(pairs < cz_pairs)
  {
    unsigned control = StabilizerSimulator::random_uint() % n_qubits;
    unsigned target = StabilizerSimulator::random_uint() % n_qubits;
    if(last_layer.size()>0)
    {      
      if (last_layer[control]== GateEnum::CZ || last_layer[target] == GateEnum::CZ)
        {
          continue;
        }
    }
    if(layer[control] == GateEnum::I && layer[target] == GateEnum::I)
    {
      layer[control] = GateEnum::CZ;
      layer[target] = GateEnum::CZ;
      pairs++;
      gates.push_back(std::make_shared<Gate>(control, target));
    }
  }
  return layer;
}

void add_single_qubit_gates(unsigned n_qubits, layer_t& last_layer, layer_t& this_layer, std::vector<bool>& has_a_t, std::vector<gate_ptr>& gates, bool old_style)
{
  for(unsigned i=0; i<n_qubits; i++)
  {
    if(this_layer[i] == GateEnum::CZ)
    {
      continue;
    }
    if(old_style)
    {
      if(!has_a_t[i])
      {
        has_a_t[i] = true;
        this_layer[i] = GateEnum::Rot;
        gates.push_back(std::make_shared<Rotation>(i));
      }
      else if(last_layer[i] == GateEnum::CZ)
      {
        GateEnum gtype;
        double r = StabilizerSimulator::random_double();
        if(r<1./3.)
        {
          gtype = GateEnum::SX;
        }
        else if(r<2./3.)
        {
          gtype = GateEnum::Rot;
        }
        else
        {
          gtype = GateEnum::SY;
        }
        this_layer[i] = gtype;
        if(gtype == GateEnum::Rot)
        {
          gates.push_back(std::make_shared<Rotation>(i));
        }
        else
        {
          gates.push_back(std::make_shared<Gate>(i, gtype));
        }
      }
    }
    else
    {
      if(last_layer[i] == GateEnum::SX || last_layer[i] == GateEnum::SY)
      {
        this_layer[i] = GateEnum::Rot;
        gates.push_back(std::make_shared<Rotation>(i));
      }
      else if(last_layer[i] == GateEnum::CZ)
      {
        GateEnum gtype;
        if(StabilizerSimulator::random_bit())
        {
          gtype = GateEnum::SX;
        }
        else
        {
          gtype = GateEnum::SY;
        }
        this_layer[i] = gtype;
        gates.push_back(std::make_shared<Gate>(i, gtype));
      }
    }
  }
}

std::vector<gate_ptr> google_circuit(unsigned n_qubits, unsigned n_layers, std::vector<std::vector<cz_t>>& cz_schema, bool old_style=false, int seed=-1)
{

  if (seed < 0)
  {
    unsigned rng_seed = std::time(nullptr);
    StabilizerSimulator::init_rng(rng_seed, 0);
  }
  else
  {
    StabilizerSimulator::init_rng(seed, 0);
  }
  std::vector<gate_ptr> flattened_circuit;
  std::vector<bool> has_a_t(n_qubits, false);
  for(unsigned i=0; i<n_qubits; i++)
  {
    flattened_circuit.push_back(std::make_shared<Gate>(i, GateEnum::H));
  }
  layer_t last_layer = get_cz_layer(n_qubits, 0, cz_schema, flattened_circuit);
  for(unsigned i=0; i<n_qubits; i++)
  {
    if(last_layer[i] == GateEnum::CZ)
    {
      continue;
    }
    if(old_style)
    {
      last_layer[i] = GateEnum::Rot;
      has_a_t[i] = true;
      flattened_circuit.push_back(std::make_shared<Rotation>(i));
    }
  }
  for(unsigned i=1; i<n_layers; i++)
  {
    unsigned layer_index = i%cz_schema.size();
    layer_t this_layer = get_cz_layer(n_qubits, layer_index, cz_schema, flattened_circuit);
    add_single_qubit_gates(n_qubits, last_layer, this_layer, has_a_t, flattened_circuit, old_style);
    last_layer = this_layer;
  }
  for(unsigned i=0; i<n_qubits; i++)
  {
    flattened_circuit.push_back(std::make_shared<Gate>(i, GateEnum::H));
  }
  return flattened_circuit;
}

std::vector<gate_ptr> infd_google_circuit(unsigned n_qubits, unsigned n_layers, double cz_fraction, bool old_style=false, int seed=-1)
{
  if (seed < 0)
  {
    unsigned rng_seed = std::time(nullptr);
    StabilizerSimulator::init_rng(rng_seed, 0);
  }
  else
  {
    StabilizerSimulator::init_rng(seed, 0);
  }
  std::vector<gate_ptr> flattened_circuit;
  std::vector<bool> has_a_t(n_qubits, false);
  for(unsigned i=0; i<n_qubits; i++)
  {
    flattened_circuit.push_back(std::make_shared<Gate>(i, GateEnum::H));
  }
  layer_t last_layer = get_cz_layer(n_qubits, cz_fraction, flattened_circuit);
  for(unsigned i=0; i<n_qubits; i++)
  {
    if(last_layer[i] == GateEnum::CZ)
    {
      continue;
    }
    if(old_style)
    {
      last_layer[i] = GateEnum::Rot;
      has_a_t[i] = true;
      flattened_circuit.push_back(std::make_shared<Rotation>(i));
    }
  }
  for(unsigned i=1; i<n_layers; i++)
  {
    layer_t this_layer = get_cz_layer(n_qubits, cz_fraction, flattened_circuit, last_layer);
    add_single_qubit_gates(n_qubits, last_layer, this_layer, has_a_t, flattened_circuit, old_style);
    last_layer = this_layer;
  }
  for(unsigned i=0; i<n_qubits; i++)
  {
    flattened_circuit.push_back(std::make_shared<Gate>(i, GateEnum::H));
  }
  return flattened_circuit;
}

// Methods for loading CZ schema from JSON

void from_json(const json_t& js, std::vector<cz_t>& cz_layer)
{
  for(auto& pair : js)
  {
    cz_layer.push_back({pair[0], pair[1]});
  }
}

json_t qasm_instructions(const std::vector<gate_ptr>& circuit, unsigned n_qubits)
{
  json_t js;
  std::vector<unsigned> all_qubits;
  for(unsigned i=0; i<n_qubits; i++)
  {
    all_qubits.push_back(i);
  }
  for(auto gate: circuit)
  {
    switch(gate->type)
    {
      case GateEnum::Rot:
      {
        json_t instruction;
        std::shared_ptr<Rotation> r = std::dynamic_pointer_cast<Rotation>(gate);
        if(r->rotation_operator.weight() > 1)
        {
          throw std::runtime_error("Unable to serialise Clifford-Recompiled circuits.");
        }
        instruction["name"] = "u1";
        instruction["params"] = {r->rotation_angle};
        instruction["texparams"] = {std::to_string(r->rotation_angle)};
        instruction["qubits"] = {r->qubit};
        instruction["memory"] = {};
        js.push_back(instruction);
        break;
      }
      case GateEnum::H:
      {
        json_t instruction;
        instruction["name"] = "h";
        instruction["params"] = {};
        instruction["texparams"] = {};
        instruction["qubits"] = {gate->qubit};
        instruction["memory"] = {};
        js.push_back(instruction);
        break;
      }
      case GateEnum::SX:
      {
        json_t instruction_1;
        instruction_1["name"] = "h";
        instruction_1["params"] = {};
        instruction_1["texparams"] = {};
        instruction_1["qubits"] = {gate->qubit};
        instruction_1["memory"] = {};
        js.push_back(instruction_1);
        json_t instruction_2;
        instruction_2["name"] = "s";
        instruction_2["params"] = {};
        instruction_2["texparams"] = {};
        instruction_2["qubits"] = {gate->qubit};
        instruction_2["memory"] = {};
        js.push_back(instruction_2);
        json_t instruction_3;
        instruction_3["name"] = "h";
        instruction_3["params"] = {};
        instruction_3["texparams"] = {};
        instruction_3["qubits"] = {gate->qubit};
        instruction_3["memory"] = {};
        js.push_back(instruction_3);
        break;
      }
      case GateEnum::SY:
      {
        json_t instruction_1;
        instruction_1["name"] = "s";
        instruction_1["params"] = {};
        instruction_1["texparams"] = {};
        instruction_1["qubits"] = {gate->qubit};
        instruction_1["memory"] = {};
        js.push_back(instruction_1);
        json_t instruction_2;
        instruction_2["name"] = "h";
        instruction_2["params"] = {};
        instruction_2["texparams"] = {};
        instruction_2["qubits"] = {gate->qubit};
        instruction_2["memory"] = {};
        js.push_back(instruction_2);
        json_t instruction_3;
        instruction_3["name"] = "s";
        instruction_3["params"] = {};
        instruction_3["texparams"] = {};
        instruction_3["qubits"] = {gate->qubit};
        instruction_3["memory"] = {};
        js.push_back(instruction_3);
        json_t instruction_4;
        instruction_4["name"] = "h";
        instruction_4["params"] = {};
        instruction_4["texparams"] = {};
        instruction_4["qubits"] = {gate->qubit};
        instruction_4["memory"] = {};
        js.push_back(instruction_4);
        json_t instruction_5;
        instruction_5["name"] = "sdg";
        instruction_5["params"] = {};
        instruction_5["texparams"] = {};
        instruction_5["qubits"] = {gate->qubit};
        instruction_5["memory"] = {};
        js.push_back(instruction_5);
        break;
      }
      case GateEnum::CZ:
      {
        json_t instruction;
        instruction["name"] = "cz";
        instruction["params"] = {};
        instruction["texparams"] = {};
        instruction["qubits"] = {gate->qubit, gate->target};
        instruction["memory"] = {};
        js.push_back(instruction);
        break;
      }
      default:
      {
        break;
      }
    }
    json_t instruction;
    instruction["name"] = "barrier";
    instruction["params"] = {};
    instruction["texparams"] = {};
    instruction["qubits"] = all_qubits;
    instruction["memory"] = {};
  }
  for(unsigned i=0; i<n_qubits; i++)
  {
    json_t instruction;
    instruction["name"] = "measure";
    instruction["params"] = {};
    instruction["texparams"] = {};
    instruction["qubits"] = {i};
    instruction["memory"] = {i};
    js.push_back(instruction);
  }
  return js;
}

json_t load_json_from_file(std::string file_name)
{
    json_t js;
    std::ifstream ifile;
    ifile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try
    {
        ifile.open(file_name);
    }
    catch (std::exception &e)
    {
      throw std::runtime_error(std::string("no such file or directory"));
    }
    ifile >> js;
    ifile.close();
    return js;
}

void dump_to_qobj(std::string experiment_name, unsigned n_qubits, std::vector<gate_ptr>& gates)
{
  json_t skeleton_qobj = load_json_from_file("./skeleton_qobj.json");
  skeleton_qobj["config"]["n_qubits"] = n_qubits;
  skeleton_qobj["config"]["memory_slots"] = n_qubits;
  json_t instruction_array = qasm_instructions(gates, n_qubits);
  skeleton_qobj["experiments"][0]["instructions"] = instruction_array;
  for(unsigned i=0; i<n_qubits; i++)
  {
    skeleton_qobj["experiments"][0]["header"]["qubit_labels"].push_back({"q0", i});
    skeleton_qobj["experiments"][0]["header"]["clbit_labels"].push_back({"c0", i});
  }
  skeleton_qobj["experiments"][0]["header"]["qreg_sizes"].push_back({"q0", n_qubits});
  skeleton_qobj["experiments"][0]["header"]["creg_sizes"].push_back({"c0", n_qubits});
  skeleton_qobj["experiments"][0]["header"]["name"] = experiment_name;
  std::string qasm_stub = skeleton_qobj["experiments"][0]["header"]["compiled_circuit_qasm"];
  for(auto gp : gates)
  {
    qasm_stub += gp->as_qasm();
    qasm_stub += "barrier ";
    for(unsigned i=0; i<n_qubits; i++)
    {
      qasm_stub += "q0["+std::to_string(i)+"]";
      if(i != n_qubits - 1)
      {
        qasm_stub += ",";
      }
    }
    qasm_stub+=";\n";
  }
  for(unsigned i=0; i<n_qubits; i++)
  {
    std::string index = std::to_string(i);
    qasm_stub += "measure q0[" + index + "] -> c0[" + index + "];\n"; 
  }
  skeleton_qobj["experiments"][0]["header"]["compiled_circuit_qasm"] = qasm_stub;
  skeleton_qobj["experiments"][0]["config"]["memory_slots"] = n_qubits;
  skeleton_qobj["experiments"][0]["config"]["n_qubits"] = n_qubits;
  std::string ofile_name = experiment_name+"_qobj.json";
  std::ofstream ofile;
  ofile.open(ofile_name);
  ofile << skeleton_qobj.dump(4) << std::endl;
  ofile.close();
}