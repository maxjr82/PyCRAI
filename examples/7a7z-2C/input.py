import sys
sys.path.append('../../')
import numpy as np
import color_rings as cr
import glob
import matplotlib.pyplot as plt

strings = np.loadtxt('energies.dat', delimiter=',', usecols=(1,7), dtype=np.str)
values = np.loadtxt('energies.dat', delimiter=',', usecols=(3,6), dtype=np.float64)
#filelist = glob.glob(r'geom*.xyz')

for i in range(len(strings)):
  xyzfile = "mol_" + strings[i][1] + "_" + strings[i][0] + ".xyz"

  if values[i][0] == 0.000:
    info = "GS, E = " + str(values[i][0]) + " eV, d$_{CC}$ = " + str(values[i][1]) + " $\AA$"
  else:
    info = "ES, E = " + str(values[i][0]) + " eV, d$_{CC}$ = " + str(values[i][1]) + " $\AA$"

  ring_size = [5,6,8]
  coords,labels,connec = cr.ring_connec(xyzfile,ring_size)

  print(len(connec))

  fname = "frm" + str(i) + "_" + xyzfile.split("_")[1] + "_" + (xyzfile.split("_")[2]).split(".")[0] 
  homa = cr.calc_homa(connec,coords,labels)

  fig, ax = cr.draw_2D_structure(coords,connec,homa,save_fig=False,num_ring=False)
  plt.tight_layout()
  cbar = fig.get_axes()[-1]
  cbar.remove()
  fig.text(0.43, 0.03, info, ha='center', fontsize=14)
  fig.tight_layout()
  fig.savefig(fname+".png",bbox_inches="tight",dpi=400)

  print("Done with {}".format(xyzfile))
  print("-------")

  plt.close()
