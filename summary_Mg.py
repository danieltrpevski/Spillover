import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import seaborn as sns

sns.set(font_scale = 1.5)
sns.set_style("ticks")
# fitting function
def eq_5_fitting_mM(V,conc,alpha, eta):
    return 1/(1+conc*eta*np.exp(alpha*V)) 


simulated_voltage=np.linspace(-250,150, 1000)
simulated_concentrations=np.logspace(-5, 1, num=35)

plt.figure(figsize=(7,8))
plt.plot(simulated_concentrations, 1/-0.0413 * np.log(1/(1.327*simulated_concentrations)), '-', color='#4daf4a', label='0.5 mM - Nowak et al., 1984', linewidth=3)
plt.plot(simulated_concentrations, 1/-0.062 * np.log(1/(0.28*simulated_concentrations)), '-', color='#e41a1c', label='1 mM - Jahr and Stevens, 1990', linewidth=3)
plt.plot(simulated_concentrations,1/-0.0486 * np.log(1/(0.49*simulated_concentrations)), '-' , color='#ff7f00', label='30 $\mu$M - Chen et al., 1992', linewidth=3)
plt.plot(simulated_concentrations,1/-0.058 * np.log(1/(0.418*simulated_concentrations)), '-' , color='#984ea3', label='2 mM - McMenimen et al., 2006', linewidth=3)
plt.plot(simulated_concentrations,1/-0.046 * np.log(1/(3.3*simulated_concentrations)), '-' , color='#377eb8', label='0.2 mM - McMenimen et al., 2006', linewidth=3)
plt.xscale('log')
plt.legend()
#plt.title('Relationship between' + '\n' + 'the voltage at which the conductance is half mazimal' + '\n' + 'and magnesium concentration', fontsize = '14')
plt.xlabel('Mg$^{2+}$ [mM]')
plt.ylabel('Voltage at g(V) = 0.5 [mV]')
#plt.grid(color='grey', linestyle='-', linewidth='0.1')
plt.savefig('v50_Mg_log.pdf')

plt.figure(figsize=(7,8))
plt.plot(simulated_concentrations, 1/-0.0413 * np.log(1/(1.327*simulated_concentrations)), '-', color='#4daf4a', label='0.5 mM - Nowak et al., 1984', linewidth=3)
plt.plot(simulated_concentrations, 1/-0.062 * np.log(1/(0.28*simulated_concentrations)), '-', color='#e41a1c', label='1 mM - Jahr and Stevens, 1990', linewidth=3)
plt.plot(simulated_concentrations,1/-0.0486 * np.log(1/(0.49*simulated_concentrations)), '-' , color='#ff7f00', label='30 $\mu$M - Chen et al., 1992', linewidth=3)
plt.plot(simulated_concentrations,1/-0.058 * np.log(1/(0.418*simulated_concentrations)), '-' , color='#984ea3', label='2 mM - McMenimen et al., 2006', linewidth=3)
plt.plot(simulated_concentrations,1/-0.046 * np.log(1/(3.3*simulated_concentrations)), '-' , color='#377eb8', label='0.2 mM - McMenimen et al., 2006', linewidth=3)
plt.legend()
#plt.title('Relationship between' + '\n' + ' the voltage at which the conductance is half mazimal' + '\n' + 'and magnesium concentration')
plt.xlabel('Mg$^{2+}$ [mM]')
plt.ylabel('Voltage at g(V) = 0.5 [mV]')
#plt.grid(color='grey', linestyle='-', linewidth='0.1')
plt.savefig('v50_Mg.pdf')

plt.figure(figsize=(7,8))
plt.plot(simulated_voltage, eq_5_fitting_mM(simulated_voltage,0.5,-0.0413, 1.327), '-', color='#4daf4a', label='0.5 mM - Nowak et al., 1984', linewidth=3)
plt.plot(simulated_voltage, eq_5_fitting_mM(simulated_voltage,1,-0.062, 0.28), '-', color='#e41a1c', label ='1 mM - Jahr and Stevens, 1990', linewidth=3)
plt.plot(simulated_voltage, eq_5_fitting_mM(simulated_voltage,0.03,-0.0486, 0.49), '-' , color='#ff7f00', label='30 $\mu$M - Chen et al., 1992', linewidth=3)
plt.plot(simulated_voltage, eq_5_fitting_mM(simulated_voltage,2,-0.058, 0.418) , '-' , color='#984ea3', label='2 mM - McMenimen et al., 2006', linewidth=3)
plt.plot(simulated_voltage, eq_5_fitting_mM(simulated_voltage,0.2,-0.046, 3.3), '-' , color='#377eb8', label='0.2 mM - McMenimen et al., 2006', linewidth=3)
plt.legend()
#plt.title('Conductance-Voltage relationship', fontsize = '14')
plt.xlabel('Voltage [mV]')
plt.ylabel('g(V)')
plt.axhline(y = 0.5, xmin=-250, xmax=150, color='grey', linestyle='-', linewidth='0.1')
#plt.grid(color='grey', linestyle='-', linewidth='0.1')
plt.savefig('g_V.pdf')
plt.show()

#import pdb
#pdb.set_trace()

