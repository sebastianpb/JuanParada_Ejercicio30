import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

x=np.loadtxt("adveccion_0.dat")
print(x[:,1])

t_max = [0,0.5,1,2]
fname=[]
for i in range(len(t_max)):
    fname.append('adveccion_%d.dat'%(t_max[i]*10))
    sol=np.loadtxt(fname[i])
    plt.plot(sol[:,1],sol[:,0])
plt.savefig('adveccion.pdf')
