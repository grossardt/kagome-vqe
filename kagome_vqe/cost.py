"""Cost function."""

def error_mitigated_cost(
    estimator, hamiltonian, ansaetze, rots, cnot_mults, theta1, theta2 = None, timeout = 3600
):
    """Cost function.
    Determines the energy of the ansatz for the given parameters using circuit estimation.
    """
    nr = len(rots)
    nm = len(cnot_mults)
    n = nr*nm
    if nm == 1:
        cnot_mitigation = False
    elif nm == 3:
        cnot_mitigation = True
    else:
        raise Exception('CNOT error mitigation is currently only supported with multiplicators 1,3,5')
    if theta2:
        runs = 2
        circuits = ansaetze * 2
        parameter_values = ([theta1]*n) + ([theta2]*n)
        observables = [hamiltonian]*(2*n)
    else:
        runs = 1
        circuits = ansaetze
        parameter_values = [theta1]*n
        observables = [hamiltonian]*n
    job = estimator.run(circuits=circuits, parameter_values=parameter_values, observables=observables)
    res = job.result()
    ev = [0]*runs
    for run in range(runs):
        evs = res.values[run*n:(run+1)*n]
        if cnot_mitigation:
            for i in rots:
                # determine the value x0 of a quadratic function y = x0 + a x + b x^2
                # through points (1,y1), (3,y3), (5,y5)
                ev[run] += (15*evs[i]-10*evs[i+1]+3*evs[i+2])/8
        else:
            ev[run] = sum(evs)
        ev[run] /= nr
    if runs == 2:
        return (ev[0],ev[1])
    return ev[0]
