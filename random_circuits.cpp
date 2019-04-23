#include "google_circuits.hpp"

#include <pybind11/pybind11.h>
#include <pybind11.stl.h>

using namespace py = pybind11;

using sqg_t = std::tuple<std::string, unsigned>;
using tqg_t = std::tuple<std::string, unsigned, unsigned>;

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
}

PYBIND11_MODULE(random_circuits, mod)
{

}