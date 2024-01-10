# a demo for computing localized homology

from blowup_complex import ProductSimplex, BlowupComplex
from complex import DelaunayComplex
from simplicial_cover import HyperCube,HyperCubeTri
import numpy as np
from matplotlib import pyplot as plt

import util

if __name__=='__main__':
    """
    linear chain with two covers 

    ```
    X:    0-1-2-3
    X0:   0-1-2  
    X1:     1-2-3
    ```
    The covers are glued by a 1D standard simplex $\Delta^1$ 
    ```
    D[1]:   0-1
    ```
    """
    
    X0 = [[0],[1],[2],[0,1],[1,2]]
    X1 = [[1],[2],[3],[1,2],[2,3]]   

    blowup = BlowupComplex([X0,X1])
    blowup.compute_persistence()
    
    
   
    """ 
    a cycle with three covers 
    ```
     X:    0--1--2
           |     |
           5--4--3
     X0: (0,1,2,3)
     X1: (2,3,4,5)
     X2: (0,1,4,5)
     ```
    The covers are glued by a 2D standard simplex $\Delta^2$ 
    """
    
    X0= [[0],[1],[2],[3],[0,1],[1,2],[2,3]]  
    X1=[[2],[3],[4],[5] ,[2,3],[3,4],[4,5]]   
    X2= [[0],[1],[4],[5],[0,1],[0,5],[4,5]]
    blowup = BlowupComplex([X0,X1,X2])
    blowup.compute_persistence(verbose=False )
    print('all basis',blowup.cycle_basis)
    

    """
    point cloud with hypercube cover including face
    """
    circles_data = util.threecircles(N=25,s=42)
    delaunay_complex = DelaunayComplex(np.array(circles_data ),4)#生成2Dmesh
    
    N=3 # columns
    M=3 # rows
    hypercube = HyperCubeTri(circles_data,delaunay_complex.edge_index(), 
                N,  M,  r=0.3)
    print("Create hyperCube cover of size:", len(hypercube.cover))
    hypercube.draw_cover()

    print(hypercube.cover)
    
    # Create the figure and axes objects
    fig, axs = plt.subplots(M, N, figsize=(N * 4, M * 3))

    # Loop over the rows and columns
    for i in range(M): #row
        for j in range(N): #col
            # Determine the current axes
            if M > 1 and N>1:
                ax = axs[ i,j]
            elif N>1:
                ax = axs[j] 
            elif M>1:
                ax = axs[i]
              
            hypercube.draw_cover_subset(i,j,ax)


    plt.tight_layout()
    plt.show()
    
    """
    compute localized homology and visualize basis
    """

    # to compute the localized homology, directly use the cover as the input of BlowupComplex()
    blowup = BlowupComplex(hypercube.cover)
    blowup.compute_persistence(verbose=False,show_diag=False)

    # visualize the last three cycles with birth 1
    cycle_edges = blowup.get_cycle_edges_by_birth(birth=1 )
    
    print("number of cylces tracked",len(cycle_edges))

    #print( "cycles with birth 1:", cycle_edges)
    delaunay_complex.draw_chains(cycle_edges)
