#include "google_circuits.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using sqg_t = std::tuple<std::string, unsigned>;
using tqg_t = std::tuple<std::string, unsigned, unsigned>;

void circ_to_python(std::vector<gate_ptr>& circ, py::list& return_list, unsigned n_qubits)
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
        return_list.append(gstr);
    }
}

py::list random_circuit(unsigned n_qubits, unsigned n_layers, std::vector<std::vector<cz_t>>& cz_schema, bool old_style, int seed)
{
    py::list gates;
    std::vector<std::vector<cz_t>> cpp_cz_schema;
    for(auto layer : cz_schema)
    {
        std::vector<cz_t> cz_layer(layer);
        cpp_cz_schema.push_back(cz_layer);
    }
    std::vector<gate_ptr> circ = google_circuit(n_qubits, n_layers, cpp_cz_schema, old_style, seed);
    circ_to_python(circ, gates, n_qubits);
    return gates;
}

py::list random_circuit(unsigned n_qubits, unsigned n_layers, double cz_fraction, bool old_style, int seed)
{
    py::list gates;
    std::vector<gate_ptr> circ = infd_google_circuit(n_qubits, n_layers, cz_fraction, old_style, seed);
    circ_to_python(circ, gates, n_qubits);
    return gates;
}

py::tuple recompiled_circuit(unsigned n_qubits, unsigned n_layers, std::vector<std::vector<cz_t>>& cz_schema, bool old_style, int seed)
{
    py::list gates, compiled_gates;
    std::vector<std::vector<cz_t>> cpp_cz_schema;
    for(auto layer : cz_schema)
    {
        std::vector<cz_t> cz_layer(layer);
        cpp_cz_schema.push_back(cz_layer);
    }
    std::vector<gate_ptr> circ = google_circuit(n_qubits, n_layers, cpp_cz_schema, old_style, seed);
    circ_to_python(circ, gates, n_qubits);
    clifford_recompilation(circ);
    circ_to_python(circ, compiled_gates, n_qubits);
    return py::make_tuple(gates, compiled_gates);
}

py::tuple recompiled_circuit(unsigned n_qubits, unsigned n_layers, double cz_fraction, bool old_style, int seed)
{
    py::list gates, compiled_gates;
    std::vector<gate_ptr> circ = infd_google_circuit(n_qubits, n_layers, cz_fraction, old_style, seed);
    circ_to_python(circ, gates, n_qubits);
    clifford_recompilation(circ);
    circ_to_python(circ, compiled_gates, n_qubits);
    return py::make_tuple(gates, compiled_gates);
}

PYBIND11_MODULE(random_circuits, mod)
{
    mod.def("generate_circuit", 
        py::overload_cast<unsigned, unsigned, std::vector<std::vector<cz_t>>&, bool, int>(&random_circuit),
        py::arg("n_qubits"),
        py::arg("n_layers"),
        py::arg("cz_schema"),
        py::arg("old_style")=false,
        py::arg("seed")=-1
        );

    mod.def("generate_connected_circuit", 
        py::overload_cast<unsigned, unsigned, double, bool, int>(&random_circuit),
        py::arg("n_qubits"),
        py::arg("n_layers"),
        py::arg("cz_fraction"),
        py::arg("old_style")=false,
        py::arg("seed")=-1
        );


    mod.def("generate_compiled_circuit",
        py::overload_cast<unsigned, unsigned, std::vector<std::vector<cz_t>>&, bool, int>(&recompiled_circuit),
        py::arg("n_qubits"),
        py::arg("n_layers"),
        py::arg("cz_properties"),
        py::arg("old_style")=false,
        py::arg("seed")=-1
        );

    mod.def("generate_connected_compiled_circuit",
    py::overload_cast<unsigned, unsigned, double, bool, int>(&recompiled_circuit),
    py::arg("n_qubits"),
    py::arg("n_layers"),
    py::arg("cz_properties"),
    py::arg("old_style")=false,
    py::arg("seed")=-1
    );

}