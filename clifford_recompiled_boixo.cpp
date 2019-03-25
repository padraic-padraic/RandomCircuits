#include "ps_stabilizer/core.hpp"
#include "ps_stabilizer/chstabilizer.hpp"
#include "ps_stabilizer/runner.hpp"

#include <unordered_map>

enum class GateEnum {H, CZ, SX, SY, Rot};

using uint_t = StabilizerSimulator::uint_t;
using pauli_t = StabilizerSimulator::pauli_t;
using layer_t = std::unordered_map<unsigned, Gate>;

uint_t zer = 0ULL;
uint_t one = 1ULL;

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
};

struct Rotation : public Gate 
{
    pauli_t rotation_operator;
    double angle;

    Rotation() : type(GateEnum::Rot) {};
    Rotation(unsigned qubit_, double angle_) : Gate() 
    {
        qubit = qubit_;
        angle = angle_;
        rotation_operator.X = zer;
        rotation_operator.e = 0;
        rotation_operator.Z = (one << qubit_);
    };

    ~Rotation() = default;

    void commute_past(Gate& other)
    {
        switch(other.type)
        {
            case GateEnum::H
            {
                uint_t shift = (one << other.qubit);
                uint_t x_temp = (rotation_operator.X & (~shift)) ^ (rotation_operator.Z & shift);
                uint_t z_temp = (rotation_operator.Z & (~shift)) ^ (rotation_operator.X & shift);
                if(!!(rotation_operator.X & rotation_operator.Z & shfit))
                {
                    rotation_operator.e ^= 2U;
                }
                break;
            }
            case GateEnum::SX:
            {
                uint_t shift = (one << i);
                if(!!(rotation_operator.Z & shift)) //Only update Z/Y
                {
                    rotation_operator.e = (rotation_operator.e +1)%4;
                    rotation_operator.X ^= shift;
                }
                break;
            }
            case GateEnum::SY:
            {
                if(!!((rotation_operator.X ^ rotation_operator.Z) & shift)) //Only update if X or Z
                {
                    rotation_operator.e ^= (!!(rotation_operator.X & shift))*2U;
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

                }
                break;
            }
            default:
                throw std::runtime_error("Commutation should stop rather than reach this point...");
        }
    }
};

