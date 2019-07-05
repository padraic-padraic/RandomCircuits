import argparse
import datetime
import json

from circuit_generation.random_circuits import generate_circuit
from circuit_generation.random_circuits import generate_connected_circuit

from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister
from qiskit.compiler import assemble
from qiskit.providers.aer import QasmSimulator

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("backend", type=str)
    parser.add_argument("schema", type=str)
    parser.add_argument("depth", type=int)
    parser.add_argument("--cz_fraction", type=float, default=0.25)
    parser.add_argument("--save_circ", action="store_true")
    parser.add_argument("--timeout")
    args = parser.parse_args()
    dims = [int(q) for q in args.schema.split('x')]
    dimension = len(dims)
    n_qubits = 1
    for dim in dims:
        n_qubits *= dim
    if args.schema == 'indf':
        circ = generate_connected_circuit(
            n_qubits,
            args.depth,
            args.cz_fraction,
            False)
    else:
        old_style = dimension == 1
        with open('CZ_schemata/{}.json'.format(args.schema), 'r') as _f:
            schema = json.load(_f)
            layers = []
            for i in range(1, len(schema)+1):
                layers.append(schema['layer_{}'.format(i)])
        circ = generate_circuit(
            n_qubits,
            args.depth,
            layers,
            old_style)
    dstr = datetime.datetime.now().isoformat()
    if args.save_circ:
        with open('./circ_{}-{}-{}.txt'.format(args.schema,
                                               args.depth,
                                               dstr), 'w') as _f:
            for gate in circ:
                _f.write("{}\n".format(gate))
    creg = ClassicalRegister(n_qubits)
    qreg = QuantumRegister(n_qubits)
    qcirc = QuantumCircuit(qreg, creg)
    for gate in circ:
        qblock_start = gate.find('[')
        qblock_end = gate.find(']')
        if gate.startswith(('H', 'SX', 'SY')):
            qubit = int(gate[qblock_start+1:qblock_end])
            if gate.startswith('H'):
                qcirc.h(qreg[qubit])
            elif gate.startswith('SX'):
                qcirc.h(qreg[qubit])
                qcirc.s(qreg[qubit])
                qcirc.h(qreg[qubit])
            else:
                qcirc.sdg(qreg[qubit])
                qcirc.h(qreg[qubit])
                qcirc.s(qreg[qubit])
                qcirc.h(qreg[qubit])
                qcirc.s(qreg[qubit])
        elif gate.startswith('Z'):
            qubit = int(gate[qblock_start+1:qblock_end].split('@')[0])
            qcirc.t(qreg[qubit])
        else:
            control, target = [
                int(q) for q in gate[qblock_start+1:qblock_end].split(',')]
        qcirc.barrier(qreg)
    qcirc.measure(qreg, creg)
    if args.timeout is not None:
        hours, minutes, seconds = (
                int(val) for val in args.timeout.strip().split(':')
            )
        timeout = datetime.timedelta(hours=hours, minutes=minutes,
                                     seconds=seconds).total_seconds()
    else:
        timeout = 18000
    qobj = assemble(qcirc, backend=QasmSimulator(), shots=5000)
    job = QasmSimulator().run(qobj, backend_options={
            'method': args.backend,
            'extended_stabilizer_measure_sampling': True,
            'extended_stabilizer_mixing_time': 3000
        })
    ostr = '{}-{}-{}.json'.format(args.schema,
                                  args.depth,
                                  dstr)
    try:
        result = job.result(timeout=timeout)
        output = result.to_dict()
    except TimeoutError:
        print("Timeout error.")
        output = {
            'success': False,
            'reason': 'timeout'
        }
    except:
        print("General error.")
        output = {
            'success': False,
            'reason': 'Uncertain.'
        }
    with open(ostr, 'w') as _f:
        json.dump(output, _f, indent=4)
