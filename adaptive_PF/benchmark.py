from mizore.method.evolution._te_benchmarker import TimeEvolutionBenchmarker, read_project_for_benchmark


def benchmark(project_path, trotter=None, exact_state_generated=False):
    if trotter is None:
        trotter = [{"n_delta_t": 1, "n_trotter_step": 1}, {"n_delta_t": 1, "n_trotter_step": 2}]
    project = read_project_for_benchmark(project_path)
    benchmarker = TimeEvolutionBenchmarker(*project)
    if not exact_state_generated:
        benchmarker.generate_exact_states()
    for item in trotter:
        benchmarker.generate_trotter_states(item)
    benchmarker.calc_infidelity_list()


if __name__ == '__main__':
    benchmark()
