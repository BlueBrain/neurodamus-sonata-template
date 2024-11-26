from neuron import h, gui
import numpy as np

h.load_file("cell.hoc")
object = h.cell()

print(h.psection())

h.cvode.active(0)

h.dt = 0.025
h.celsius = 32.0
h.tstop = 2000
h.v_init = -70


stim = [
    h.IClamp(0.5, sec=object.soma[0]),
    h.IClamp(0.5, sec=object.soma[0]),
    h.IClamp(0.5, sec=object.soma[0]),
]

stim[0].delay = 100
stim[0].dur = 1500
stim[0].amp = 0.01  # 10pA

stim[1].delay = 1700
stim[1].dur = 1500
stim[1].amp = 0.016  # 16pA

stim[2].delay = 3300
stim[2].dur = 1500
stim[2].amp = 0.022  # 22pA

h('load_file("vm.ses")')

sec = object.soma[0]

h.define_shape()


def initialize():
    h.finitialize()
    h.run()


initialize()
