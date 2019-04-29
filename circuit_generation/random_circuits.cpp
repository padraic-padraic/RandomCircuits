#include "google_circuits.hpp"

#include <pybind11/pybind11.h>
#include <pybind11.stl.h>

using namespace py = pybind11;

using sqg_t = std::tuple<std::string, unsigned>;
using tqg_t = std::tuple<std::string, unsigned, unsigned>;

void circ_to_python(std::vector<gate_pt>& circ, py::list& return_list)
{
    for (auto gptr : circ)
    {
        std::string gstr;
        if (gptr->type == GateEnum::Rot)
        {
            std::shared_ptr<Rotation> r = std::dynamic_pointer_cast<Rotation>(gptr);
            gstr = r->to_string(n_qubits);
        }
        else
        {
            gstr = *gptr; //Can use implicit conversion here
        }
        gates.append(gstr);
    }
}

py::list random_circuit(unsigned n_qubits, unsigned n_layers, py::list& cz_schema, bool old_style)
{
    py::list gates;
    std::vector<std::vector<cz_t>> cpp_cz_schema;
    for(auto layer : cz_schema)
    {
        std::vector<cz_t> layer;
        for(auto pair: layer)
        {
            layer.push_back({pair[0], pair[1]});
        }
        cz_schema.push_back(layer);
    }
    std::vector<gate_ptr> circ = google_circuit(n_qubits, n_layers, cpp_cz_schema, old_style);
    circ_to_python(circ, gates);
    return gates;
}

py::list random_circuit(unsigned n_qubits, unsigned n_layers, double cz_fraction, bool old_style)
{
    py::list gates;
    std::vector<gate_ptr> circ = infd_google_circuit(n_qubits, n_layers, cz_fraction, old_style);
    circ_to_python(circ, gates);
    return gates;
}

py::tuple recompiled_circuit(unsigned n_qubits, unsigned n_layers, py::list& cz_schema, bool old_style)
{
    py::list gates, compiled_gates;
    std::vector<std::vector<cz_t>> cpp_cz_schema;
    for(auto layer : cz_schema)
    {
        std::vector<cz_t> layer;
        for(auto pair: layer)
        {
            layer.push_back({pair[0], pair[1]});
        }
        cz_schema.push_back(layer);
    }
    std::vector<gate_ptr> circ = google_circuit(n_qubits, n_layers, cpp_cz_schema, old_style);
    circ_to_python(circ, gates);
    clifford_recompilation(circ);
    circ_to_python(circ, compiled_gates);
    return py::make_tuple(gates, compiled_gates);
}

py::tuple recompiled_circuit(unsigned n_qubits, unsigned n_layers, double cz_fraction, bool old_style)
{
    py::list gates, compiled_gates;
    std::vector<gate_ptr> circ = infd_google_circuit(n_qubits, n_layers, cz_fraction, old_style);
    circ_to_python(circ, gates);
    clifford_recompilation(circ);
    circ_to_python(circ, compiled_gates);
    return py::make_tuple(gates, compiled_gates);
}

PYBIND11_MODULE(random_circuits, mod)
{
    mod.def("generate_circuit", &random_circuit,
        py::arg("n_qubits"),
        py::arg("n_layers"),
        py::arg("cz_properties"),
        py::arg("old_style")=false
        );
}