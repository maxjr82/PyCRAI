import matplotlib.cm as cm
import numpy as np
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import numpy.linalg as linalg
import sys

def ring_connec(xyzfile,ring_size,find_ring=True):
  """
    Function that detects rings within a 2D structure

    Args:
        xyzfile (str): name of the .xyz file that contain the molecular coordinates
        ring_size ([nx1] array or list): array/list containing the size of the rings to be searched
        find_ring (bool): find ring or just read coordinates/labels
    Returns:
        connec ([nrings,ring_size] array): array containing the connectivity of all rings
        coords/labels: labels and coordinates of the atoms

  __authors__ = "Leonardo A. Cunha, Max Pinheiro Jr."
  __date__ = "2018-07-11"
  """

  patches = []
  # Reading coordinates and labels of all atoms in the .xyz file
  coords = np.genfromtxt(xyzfile,skip_header=2,usecols=(1,2,3))
  labels = np.genfromtxt(xyzfile, unpack=True, usecols=0, skip_header=2, dtype=str)
  if (not find_ring):
    return coords,labels
  # Building the Graph that contains the connectivities of all atoms, except hydrogen
  graph = {}
  for atom1,name in enumerate(labels):
    if name != 'H':
      graph[atom1] = []
      for atom2 in range(len(coords)):
        if (atom1 != atom2 and linalg.norm(coords[atom1]-coords[atom2])>1.30 and linalg.norm(coords[atom1]-coords[atom2])<1.70):
            graph[atom1].append(atom2)

  # Use DFS to search for the rings 
  cycles = [[node]+path for node in graph for path in dfs(graph, node, node,ring_size)]  

  # Eliminate possible duplicates and unwanted rings
  temp = []
  for i in cycles:
   if (len(i)-1) in ring_size:
       temp.append(i)
  hexa = [temp[0]]
  for temp1 in (temp[1:]):
      for temp2 in temp:
        hexa_comp = [set(i) for i in hexa]
        if set(temp1)-set(temp2) != set([]):
               if set(temp2) not in hexa_comp:
                   hexa.append(temp2)
            
  connec = np.asarray(hexa)
  return coords,labels,connec

def dfs(graph, start, end,ring_size):
  """
   Depth First Search (DFS) algorithm to search paths in a graph
   
   Args:
        graph (dictionary of lists)
        start,end (int)
        ring_size (list or array)

    Returns:
        path
  """

  fringe = [(start, [])]
  while fringe:
      state, path = fringe.pop()
      if len(path)>(max(ring_size)+1):
        continue
      if len(path) in ring_size and state == end:
          yield path
          continue
      for next_state in graph[state]:
          if next_state in path:
              continue
          fringe.append((next_state, path+[next_state]))

def calc_homa(rings,coords,labels):
  """
    Calculates the 'Harmonic Oscillator Model of Aromaticiy' (HOMA) inde

    Args:
        rings (array): the atom connectivity of the rings
        coords ([natom,3] array): atomic cartesian coordinates 

    Returns:
        homa (array): the HOMA index for each ring

  __authors__ = "Leonrdo A. Cunha, Max Pinheiro Jr."

  """

  homa = []
  Ropt_CC = 1.388
  Ropt_CN = 1.334
  Ropt_NN = 1.309
  alpha_CC = 257.7
  alpha_CN = 93.52
  alpha_NN = 130.33

  for n in rings:
    paircontr = 0
    R = 0
    for i in range(len(n)-1):
      atom=n[i]
      atom_next=n[(i+1)]
      R = linalg.norm(coords[atom]-coords[atom_next])
      if labels[atom] == 'C' and labels[atom_next] == 'C':
        Ropt = Ropt_CC
        alpha = alpha_CC
      elif ((labels[atom] == 'C' and labels[atom_next] == 'N') or (labels[atom] == 'N' and labels[atom_next] == 'C')):
        Ropt = Ropt_CN
        alpha = alpha_CN
      elif labels[atom] == 'N' and labels[atom_next] == 'N':
        Ropt = Ropt_NN
        alpha = alpha_NN
      paircontr += (alpha/len(n))*(Ropt - R)**2
    homa.append(1-paircontr)
  return homa
 
def draw_2D_structure(coords,connec,colors=None,num_ring=False,save_fig=True,fig_name="test.png"):
  """
    Draw the 2D structure of the molecule and color each ring with a specific color

    Args:
       coords (array): atomic cartesian coordinates
       connec: ring connectivity
       colors: color for each ring
       fig_name: name of the figure
    
    Returns:
        fig,ax: matplotlib objects
  __authors__ = "Leonardo A. Cunha, Max Pinheiro Jr."
  """
  patches = []
  x_centers = []
  y_centers = []
  for ring in connec:
      coord = []
      for atom in ring:
        coord.append(list(coords[atom,:2]))
      coord = np.asarray(coord)
      xc = coord[:-1,0].sum()/(len(coord[:-1,0]))
      yc = coord[:-1,1].sum()/(len(coord[:-1,1]))
      polygon = Polygon(coord, True,ec='black',lw=20,ls='solid')
      patches.append(polygon)
      x_centers.append(xc)
      y_centers.append(yc)
  
  fig, ax = plt.subplots()
  p = PatchCollection(patches,cmap=cm.jet)
  if colors:
    p.set_array(np.array(colors))
  else:
    p.set_facecolor('none')
  p.set_edgecolor('black')
  ax.add_collection(p)
  if colors:
    cbar = fig.colorbar(p, ax=ax)
#    cbar.set_clim(vmin=-2.5,vmax=1.0)
  if num_ring:
    for i in range(1,len(connec)+1):
      ax.annotate(str(i),xy=(x_centers[i-1],y_centers[i-1]),ha="center",va="center",color="k")
  ax.autoscale_view()
  plt.axis('off')
  plt.tight_layout()
  if save_fig:
    fig.savefig(fig_name,bbox_inches='tight',dpi=300)
  else:
    return fig,ax
