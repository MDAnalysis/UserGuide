import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysisTests.datafiles import AUX_EDR
aux = mda.auxiliary.EDR.EDRReader(AUX_EDR)
temp = aux.get_data("Temperature")
plt.plot(temp["Time"], temp["Temperature"])
plt.ylabel("Temperature [K]")
plt.xlabel("Time [ps]")
plt.show()