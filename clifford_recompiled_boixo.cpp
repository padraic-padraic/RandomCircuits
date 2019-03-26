#include "ps_stabilizer/core.hpp"
#include "ps_stabilizer/chstabilizer.hpp"
#include "ps_stabilizer/runner.hpp"

#include <array>
#include <unordered_map>
#include <vector>

enum class GateEnum {H, CZ, SX, SY, Rot};

using uint_t = StabilizerSimulator::uint_t;
using pauli_t = StabilizerSimulator::pauli_t;
using complex_t = StabilizerSimulator::complex_t;
using branch_t = std::pair<double, pauli_t>;
using layer_t = std::unordered_map<unsigned, GateEnum>;
using cz_t = std::pair<unsigned, unsigned>;

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
    virtual Gate(Gate& other) : qubit(other.qubit), target(other.target), type(other.type) {};
    virtual ~Gate() = default;
    virtual void commute_past(Gate& other) { return; };
};

struct Rotation : public Gate 
{
  pauli_t rotation_operator;
  double rotation_angle;
  double branch_angle; 

  double p_threshold;
  std::arrary<complex_t, 2> coeffs;

  Rotation() : type(GateEnum::Rot), rotation_angle(M_PI/4.),
  branch_angle(M_PI/2.), p_threshold(0.5) {};
  
  Rotation(unsigned qubit_) : Rotation() 
  {
    qubit = qubit_;
    rotation_operator.X = zer;
    rotation_operator.e = 0;
    rotation_operator.Z = (one << qubit_);
  };

  ~Rotation() = default;

  void commute_past(Gate& other) override
  {
    switch(other.type)
    {
      case GateEnum::H
      {
        uint_t shift = (one << other.qubit);
        uint_t x_temp = (rotation_operator.X & (~shift)) ^ (rotation_operator.Z & shift);
        uint_t z_temp = (rotation_operator.Z & (~shift)) ^ (rotation_operator.X & shift);
        if (!!(rotation_operator.X & rotation_operator.Z & shfit))
        {
          rotation_operator.e ^= 2U;
        }
        break;
      }
      case GateEnum::SX:
      {
        uint_t shift = (one << i);
        if (!!(rotation_operator.Z & shift)) //Only update Z/Y
        {
          rotation_operator.e = (rotation_operator.e + 3)%4;
          rotation_operator.X ^= shift;
        }
        break;
      }
      case GateEnum::SY:
      {
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
        uint_t control_shift = (one << qubit);
        uint_t target_shift = (one << target);
        //Update only if we act as X or Y on one of the qubits.
        if(!!(rotation_operator.X & control_shift) || !!(rotation_operator.X & target_shift))
        {
          rotation_operator.Z ^= (!!(rotation_operator.X & control_shift)) * target_shift;
          rotation_operator.Z ^= (!!(rotation_operator.X & target_shift)) * control_shift;
          rotation_operator.e ^= (!!(rotation_operator.X & control_shift) && !!(rotation_operator.X & target_shift))*2U;
        }
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

void boixo_clifford_recompilation(std::vector<Gate>& gate_sequence)
{
  std::vector<unsigned> non_clifford_positions;
  unsigned n_gates = gate_sequence.size();
  for(unsigned i = n_gates-1; i<0; i--)
  {
    if(gate_sequence[i].type == GateEnum::Rot)
    {
      clifford_positions.push_back(i);
    }
  }
  for(auto position : non_clifford_positions)
  {
    Gate non_clifford = gate_sequence[position];
    unsigned current_pos = position;
    unsigned next_pos = position+1;
    while(next_pos < gate_sequence.size())
    {
      Gate next = gate_sequence[next_pos];
      if(next.type == GateEnum::H || next.type == GateEnum::Rot)
      {
        break;
      }
      non_clifford.commute_past(next);
      gate_sequence[current_pos] = next;
      current_pos++;
      next_pos++;
    }
    gate_sequence[current_pos] = non_clifford;
  }
}
