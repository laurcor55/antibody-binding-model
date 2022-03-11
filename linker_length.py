import matplotlib.pyplot as plt
import numpy as np

def p_r(r, n, l_p, b):
  l_c = b*n # contour length
  p1 = (3/(4*np.pi*l_p*l_c))**(3/2)
  p2 = np.exp(-3*r**2/(4*l_p*l_c))
  p3 = 5*l_p/(4*l_c) - 2*r**2/(l_c**2) + 33*r**4/(80*l_p*l_c**3) + 79*l_p**2/(160*l_c**2) + 329*r**2*l_p/(120*l_c**3) - 6799*r**4/(1600*l_c**4) + 3441*r**6/(2800*l_p*l_c**5) -1089*r**8/(12800*l_p**2*l_c**6)

  return(4 * np.pi * r**2 * p1 * p2 * (1 - p3))


n = 24 # Amino acids in linker
l_p = 4.5 # persistence length in angstroms
b = 3.8 # carbon carbon bond in angstroms
dr = 0.01 # step size to create r
r = np.arange(0, 100, dr) # contour lengths in angstroms

pdf = p_r(r, n, l_p, b)
plt.plot(r, pdf)
plt.title('Probability of Contour Length')
plt.ylabel('p(r)')
plt.xlabel('Contour Length r (Angstrom)')
plt.show()

average_contour_length = np.sum(pdf*r)*dr
print(average_contour_length)