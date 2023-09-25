from bluepysnap import Simulation
import matplotlib.pyplot as plt

sim = Simulation('run_no_dynamics/simulation_config.json')
report_soma = sim.reports['soma'].filter('All')
report_soma.trace()
plt.savefig('soma_trace.png')
print('Done')