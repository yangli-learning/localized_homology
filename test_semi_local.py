# a demo for computing localized homology

from blowup_complex import  BlowupComplex
from complex import DelaunayComplex
from simplicial_cover import HyperCube
import numpy as np
from matplotlib import pyplot as plt
from semi_blowup_complex import SemiBlowupComplex
import util
import persistence as persis
import cechmate as cm
import phat



def compute_standard_persis(ordered_simplices ): 
  ## Setup boundary matrix and reduce
  columns = cm.solver._simplices_to_sparse_pivot_column(ordered_simplices)
  print("boundary of each simplex",columns)
  boundary_matrix = phat.boundary_matrix(
      columns=columns , representation=phat.representations.sparse_pivot_column
  )
  pairs = boundary_matrix.compute_persistence_pairs()
  pairs.sort() 

  ## Setup persistence diagrams by reading off distances
  dgms = persis._process_distances(pairs,ordered_simplices,show_diag=True )
  print("persistence pairs")
  for p in pairs:
      print(p)
  print(dgms)
  dgms,pairs = persis._add_unpaired(dgms, pairs, ordered_simplices)

  print("persistence dgm")
  print(dgms)
  #print(dgms[1])
  for i in range(len(dgms)):
      persis.plotbarcode_BC(dgms,i)  


if __name__=='__main__':
  """
  test semi-local persistence homology computation
  """
  """ 
  # Example 1
  X0 = [[0],[1],[2],[1,2],[0,1],[0,2],[0,1,2] ]
  X1 = [[0],[2],[3],[0,3],[2,3],[0,2],[0,2,3]]
  dist = dict({(0,):0, (1,):0,
          (2,):1, (3,):1, (1,2):1, (0,3):1,
          (0,1):2, (2,3):2,
          (0,2):3,
          (0,2,3):4,
          (0,1,2):5})
  """
  
  # Example 2
  
  X0 = [[0],[1],[2],[0,1],[1,2] ]
  X1 = [[1],[2],[3],[1,2],[2,3] ]
  dist = dict({(0,):0, 
          (1,):1, (0,1):1,
          (2,):2, (1,2):2, 
          (3,):3, (2,3):3})
  

  # compute standard filtration based on the distance function
  compute_standard_persis( [ (list(s), t)  for s, t in dist.items()])
  

  # compute localized homology based on cover
  blowup = BlowupComplex([X0,X1] )
  blowup.compute_persistence(verbose=True,show_diag=True)

  # compute semi-localized homology: use standard filtration on dim 0 subsets of 
  # localized homology 
  semi_blowup = SemiBlowupComplex([X0,X1], dist )
  semi_blowup.compute_persistence(verbose=True,show_diag=True)
