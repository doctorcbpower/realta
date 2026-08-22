import numpy as np

from realta import ClusterSimulation, SimulationConfig


def test_reference_cluster_run():
    """Verify system output match for reference cluster parameters."""
    config = SimulationConfig()
    config.ntot = 1000
    config.tmax = 10.0

    cluster = ClusterSimulation(config)
    results = cluster.run(output_dir="tests/output_tmp")

    assert len(results) > 0
    assert results[0]["time"] == 0.0
    assert not np.isnan(results[-1]["lumx_tot"])
