# a demo for computing localized homology

from blowup_complex import ProductSimplex, BlowupComplex
from complex import DelaunayComplex
from simplicial_cover import HyperCube,HyperCubeTri
import numpy as np
from matplotlib import pyplot as plt
from datasets import *

import util


if __name__=='__main__':
    """
    Example 1: linear chain with two covers 

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
    Example 2: a cycle with three covers 
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
    Example 3: point cloud with hypercube cover including faces
    """
    circles_data = util.threecircles(N=25,s=42)
    delaunay_complex = DelaunayComplex(np.array(circles_data ),4)#生成2Dmesh
    
    N=2 # columns
    M=2 # rows
    hypercube = HyperCubeTri(circles_data,delaunay_complex.edge_index(), 
                N,  M,  r=0.3)
    print("Create hyperCube cover of size:", len(hypercube.cover))
    hypercube.draw_cover()
    hypercube.visualize_cover_subsets()
    # compute the localized homology 
    blowup = BlowupComplex(hypercube.cover)
    blowup.compute_persistence(verbose=False,show_diag=False)

    #visualize the last three cycles with birth 1
    cycle_edges = blowup.get_cycle_edges_by_birth(birth=1 )
    
    #print( "cycles with birth 1:", cycle_edges)
    delaunay_complex.draw_chains(cycle_edges)
    

    """
    Example 4: point cloud from image using hypercube cover
    """
    dataset = MNIST(root='data/MNIST', download=True, train=False)
    dataset_2d = MNIST2D(dataset,   200)
    test_point_cloud =  dataset_2d[100][0].numpy()
    complex = DelaunayComplex(np.array(test_point_cloud),2) 
    
    hypercube = HyperCubeTri(test_point_cloud, complex.edge_index(), 
              N=3, #N
              M=3, #M, the numbers of covers is N*M  
              r=0.3
              )
    hypercube.draw_cover()

    #hypercube.visualize_cover_subsets() 


    blowup = BlowupComplex(hypercube.cover)
    blowup.compute_persistence(verbose=True, show_diag=False )
    
    # visualize the last three cycles with birth 0 and 1 

    cycle_edges0 = blowup.get_cycle_edges_by_birth(birth=0 )
    complex.draw_chains(cycle_edges0)

    print( "cycles with birth 0:", cycle_edges0) 

    cycle_edges1 = blowup.get_cycle_edges_by_birth(birth=1 )
    print( "cycles with birth 1:", cycle_edges1) 
    complex.draw_chains(cycle_edges1)
 


