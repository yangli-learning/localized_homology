# a demo dor computing localized homology

from blowup_complex import ProductSimplex, BlowupComplex
from complex import DelaunayComplex
from simplicial_cover import HyperCube
import numpy as np

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
    """ 
    X0 = [[0],[1],[2],[0,1],[1,2]]
    X1 = [[1],[2],[3],[1,2],[2,3]]   

    blowup = BlowupComplex([X0,X1])
    blowup.compute_persistence()
    """

    

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
    """ 
    X0= [[0],[1],[2],[3],[0,1],[1,2],[2,3]]  
    X1=[[2],[3],[4],[5] ,[2,3],[3,4],[4,5]]   
    X2= [[0],[1],[4],[5],[0,1],[0,5],[4,5]]
    blowup = BlowupComplex([X0,X1,X2])
    blowup.compute_persistence()
    """ 

    """
    point cloud with hypercube cover
    """
    circles_data = util.threecircles(N=20,s=42)
    delaunay_complex = DelaunayComplex(np.array(circles_data ),4)#生成2Dmesh
    #print("edge list:", test_complex.edge_index()) #边

    delaunay_complex.draw_complex() #绘制
    #plt.scatter(circles_data[:,0],circles_data[:,1],color = "red")
    #plt.show() 

    a = HyperCube(circles_data,delaunay_complex.edge_index(), 
                N=6,  M=4,  r=0.8 )
    print("HyperCube cover size:", len(a.cover()))
    a.draw_cover()
    #return a list,with N*M sublists represnted N*M different covers(only vertices and edges)

    # to compute the localized homology, directly use the cover as the input of BlowupComplex()
    blowup = BlowupComplex(a.cover())
    blowup.compute_persistence(show_diag=False)