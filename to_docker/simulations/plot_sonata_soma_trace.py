from bluepysnap import Simulation
import matplotlib.pyplot as plt


checker = False
while checker == False:
    celltype = input("Enter the Cell Type (int): 0 for MC, 1 for mTC, or 3 for GC: ")
    if celltype=="0":
        sim = Simulation('run_no_dynamics/simulation_MC_config.json')
        checker = True
        Ctype = "MC"
    elif celltype=="1":
        sim = Simulation('run_no_dynamics/simulation_mTC_config.json')
        checker = True
        Ctype = "mTC"
    elif celltype=="2":
        sim = Simulation('run_no_dynamics/simulation_GC_config.json')
        checker = True
        Ctype = "GC"
    else:
        print("Please enter an integer between 0 and 2. Quotes are NOT needed.")
report_soma = sim.reports['soma'].filter('All')
report_soma.trace()
filename = Ctype+'_soma_trace.png'
plt.savefig(filename)
print('Done')
