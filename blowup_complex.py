
import util 
import itertools
import persistence as persis
import cechmate as cm  

class ProductSimplex:
    """
    class ProductSimplex:
    
    The blowup chain complex is a direct product of simplicial complex $X^J$ and $\Delta^J$
    Each simplex is a pair $(\sigma,\delta)$ where $\sigma\in X^J$ and $\delta\in\Delta^J$
    """
    def __init__(self, sigma,delta,id=0):
        self.sigma =  sigma  
        self.delta  =  delta  
        self.id = id
        self.dim_s = max(len(self.sigma)-1,0) 
        self.dim_d = max(len(self.delta)-1,0)
        self.dim = self.dim_s +  self.dim_d 
        self.boundary=[] 
        
    def print(self):
        print( "id=",self.id,"simplex=",(self.sigma,self.delta), "dim=",self.dim)
        if self.boundary:
            print("boundary:",self.boundary)
            
    def to_ordered_simplex(self):
        simplex = self.sigma if len(self.sigma)==self.dim+1 else self.sigma + [0]*(self.dim+1-len(self.sigma)) # append zeros
        return (simplex, self.dim_d)
    
class BlowupComplex():
    """
    class BlowupComplex

    Expernal Functions:
    - BlowupComplex(simplicial_cover): build blowup complex
    - compute_persistence(): compute and visualize the persistence barcode
    - print_complex(): print out information of the blowup complex
    
    Internal Functions:
    - construct_from_cover(): construct blowup complex $X^U$ as a direct product 
    - compute_boundary_delta(D_J): compute boundary of $\Delta^J$
    - compute_boundary(D_J,B_delta): compute boundary of blowup complex $X^U$

    """
    def __init__(self, simplicial_cover):
        self.U = simplicial_cover  
        self.X_U = []   # list of lists of ProductSimplex (to do : make it iterable)
        self._construct_from_cover()
        
    def print_complex(self):
        # display the constructed blowup complex X_U
        for X_U_t in self.X_U: 
            for sp in X_U_t:
                sp.print()

    
    def _construct_from_cover(self): 
        i_start=0
        D_J  = util.make_standard_simplex(len(self.U))  
        for delta in D_J:
            
            X_J = util.subset_intersection([ self.U[d] for d in delta])  
            X_U_t = self._simplicial_product(X_J,delta,i_start)
            i_start += len(X_U_t)
            self.X_U.append( X_U_t  )
        
        # map simplex tuples to global id
        flattened_list =[ ( tuple(sp.sigma),tuple(sp.delta) )  for i in range(len(D_J)) for sp in self.X_U[i]]
        self.sp_to_id  = dict( [(sp,i) for i,sp in enumerate(flattened_list)])

        B_delta  = self._compute_boundary_delta(D_J) 
        print(B_delta)
        self._compute_boundary(D_J,B_delta)
        

    def _compute_boundary_delta(self,D_J):
        ordered_simplex_delta =  [  ( delta, len(delta)-1) for delta in D_J]  
        B_delta = cm.solver._simplices_to_sparse_pivot_column(ordered_simplex_delta ,False)
        return B_delta

    def _simplicial_product(self,complex, delta,i_start):
        # e.g. complex = [[1], [2], [1, 2]]   # detla=[0,1]
        # resulting productsimplex are 
        # ProductSimplex([1],[0,1]), ProductSimplex([2],[0,1]),ProductSimplex([1,2],[0,1]) 

        return [ProductSimplex(sigma,delta,i+i_start) for i, sigma in enumerate(complex)] 

    def _compute_boundary(self,D_J,B_delta ):
        for t,delta in enumerate(D_J): 
            
            # construct a filtration on Xi
            # fill in block diagonal boundary matrix by using builtin boundary operator
            ordered_simplex = [  ( sp.sigma, sp.dim_d) for sp in self.X_U[t] ]
            B = cm.solver._simplices_to_sparse_pivot_column(ordered_simplex ,False)
            
            for i,sp in enumerate(self.X_U[t]): 
                sp.boundary = [self.X_U[t][j].id for j in B[i][1]] 
            
            # fill in the upper triangular area for deltas with non-empty boundary
            if len(delta)>1: 
                for sp in self.X_U[t]: 
                    dim,boundary_delta  = B_delta[t] 
                    
                    #[2,3]x [0,1] -> [2,3]x[0], [2,3]x[1]
                    for tt in boundary_delta: 
                        boundary_delta_simplex = D_J[tt] 
                        sp.boundary.append(  self.sp_to_id[  (tuple(sp.sigma), tuple(  boundary_delta_simplex   ))]) 
                    
            for sp in self.X_U[t]:
                sp.boundary = sorted(sp.boundary,key=lambda x: (  x  ))
                
    def compute_persistence(self,verbose=False,show_diag=True):
        # convert simplices and boundary to the column format for persistence
        ordered_simplices = [sp.to_ordered_simplex() for X_U_T in self.X_U  for sp in  X_U_T   ]  
        if verbose:
            print('simplices',ordered_simplices)
        
        columns = [ (sp.dim, sp.boundary )    for X_U_T in self.X_U for sp in  X_U_T]  

        if verbose:
            print('boundary columns')
            for col in columns:
                print(col)
        # for coloring localized cycles
        cover = [ '-'.join([str(s) for s in sp.delta ])  for X_U_T in self.X_U  for sp in  X_U_T   ]  
        if verbose:
            print('cover',cover)
        self.dgms = persis.compute_persistence(ordered_simplices, columns,cover, show_diag=show_diag) 
    