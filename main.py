# a demo dor computing localized homology

from blowup_complex import ProductSimplex, BlowupComplex

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
    blowup.compute_persistence()