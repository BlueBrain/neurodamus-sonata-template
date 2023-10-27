from bluepysnap import Simulation
import matplotlib.pyplot as plt
import sys

def check_integer(integer, arg):
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
        if arg:
            print("Wrong value, no cell associated with this value. Please enter only 0, 1, or 2")
            sim, checker, Ctype = get_arg()
        if not arg:
            print("Please enter an integer between 0 and 2. Quotes are NOT needed.")
            return("0", False, "0")
    return(sim, checker, Ctype)

def get_arg():
    checker = False
    while not checker:
        celltype = input("Enter the Cell Type (int): 0 for MC, 1 for mTC, or 2 for GC: ")
        sim, checker, Ctype = check_integer(celltype, False)
    return(sim, checker, Ctype)

if len(sys.argv)>1:
    celltype = sys.argv[1]
    sim, checker, Ctype = check_integer(celltype, True)
else:
    sim, checker, Ctype = get_arg()

report_soma = sim.reports['soma'].filter('All')
report_soma.trace()
filename = Ctype+'_soma_trace.png'
plt.savefig(filename)
print('Done')
