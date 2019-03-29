#include "ps_stabilizer/core.hpp"
#include "ps_stabilizer/chstabilizer.hpp"
#include "ps_stabilizer/rng.hpp"

#define _USE_MATH_DEFINES

#include <array>
#include <cmath>
#include <memory>
#include <unordered_map>
#include <vector>

enum class GateEnum {I, H, CZ, SX, SY, Rot};

struct Gate;

using uint_t = StabilizerSimulator::uint_t;
using pauli_t = StabilizerSimulator::pauli_t;
using complex_t = StabilizerSimulator::complex_t;
using branch_t = std::pair<double, pauli_t>;
using layer_t = std::vector<GateEnum>;
using cz_t = std::pair<unsigned, unsigned>;
using gate_ptr = std::shared_ptr<Gate>;

pauli_t IDENTITY;

uint_t zer = 0ULL;
uint_t one = 1ULL;

const complex_t pi_over_8_phase(0., M_PI/8);
const complex_t root_omega = std::exp(pi_over_8_phase);
const complex_t root_omega_star = std::conj(root_omega);

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
    rotation_operator.e = 0;
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
};

void boixo_clifford_recompilation(std::vector<gate_ptr>& gate_sequence)
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

std::vector<gate_ptr> boixo_circuit(unsigned n_qubits, unsigned n_layers, std::vector<std::vector<cz_t>>& cz_schema)
{
  unsigned seed = std::time(nullptr);
  #ifdef _OPENMP
  #pragma parallel
  {
    init_rng(seed, omp_get_thread_num());
  }
  #else
  init_rng(seed, 0);
  #endif
  std::vector<gate_ptr> flattened_circuit;
  for(unsigned i=0; i<n_qubits; i++)
  {
    flattened_circuit.push_back(std::make_shared<Gate>(i, GateEnum::H));
  }
  layer_t last_layer(n_qubits, GateEnum::I);
  for(auto cz : cz_schema[0])
  {
    last_layer[cz.first] = GateEnum::CZ;
    last_layer[cz.second] = GateEnum::CZ;
    flattened_circuit.push_back(std::make_shared<Gate>(cz.first, cz.second));
  }
  for(unsigned i=0; i<n_qubits; i++)
  {
    if(last_layer[i] == GateEnum::CZ)
    {
      continue;
    }
    last_layer[i] = GateEnum::Rot;
    flattened_circuit.push_back(std::make_shared<Rotation>(i));
  }
  for(unsigned i=1; i<n_layers; i++)
  {
    layer_t this_layer(n_qubits, GateEnum::I);
    auto czs = cz_schema[n_layers%cz_schema.size()];
    for(auto cz: czs)
    {
      this_layer[cz.first] = GateEnum::CZ;
      this_layer[cz.second] = GateEnum::CZ;
      flattened_circuit.push_back(std::make_shared<Gate>(cz.first, cz.second));
    }
    for(unsigned i=0; i<n_qubits; i++)
    {
      if(this_layer[i] == GateEnum::CZ)
      {
        continue;
      }
      if(last_layer[i] == GateEnum::SX || last_layer[i] == GateEnum::SY)
      {
        this_layer[i] = GateEnum::Rot;
        flattened_circuit.push_back(std::make_shared<Rotation>(i));
      }
      else
      {
        GateEnum gtype;
        if(random_bit())
        {
          gtype = GateEnum::SX;
        }
        else
        {
          gtype = GateEnum::SY;
        }
        this_layer[i] = gtype;
        flattened_circuit.push_back(std::make_shared<Gate>(i, gtype));
      }
    }
    last_layer = this_layer;
  }
  return flattened_circuit;
}
