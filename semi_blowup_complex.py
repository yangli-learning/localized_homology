import persistence as persis
from blowup_complex import ProductSimplex, BlowupComplex
import cechmate as cm  
import numpy as np

class MetricProductSimplex(ProductSimplex):
    """
    A subclass of ProductSimplex that include a distance parameter `dist`
    which represents the inclusion time of the sigma in the
    underlying filtered simplicial complex
    
    """
    def __init__(self,sigma,delta,dist,id=0):
        super().__init__(sigma, delta,id)
        self.dist = dist # distance parameter 

    def print(self):
        super().print()
        print('dist:', self.dist)
    
    def to_ordered_simplex(self,max_dist):
        simplex = self.sigma if len(self.sigma)==self.dim+1 else self.sigma + [0]*(self.dim+1-len(self.sigma)) # append zeros
        if self.dim_d == 0:
            t = self.dist
        else:
            t = max_dist + self.dim_d
        return (simplex, t )

class SemiBlowupComplex(BlowupComplex):
    def __init__(self, simplicial_cover,dist ):
        self.dist= dist # a dictionary containing distance for each simplex  

        super().__init__(simplicial_cover)
    
    def _simplicial_product(self, complex, delta, i_start): 
        
        return [MetricProductSimplex(sigma,delta,self.dist[tuple(sigma)],i+i_start) for i, sigma in enumerate(complex)] 
    
    def compute_persistence(self, verbose=False, show_diag=True,
                             compute_basis=True # not implemented
                             ):
        
        # modify the distance for simplices
        ordered_simplices =  [sp.to_ordered_simplex(max(self.dist.values())) for X_U_T in self.X_U  for sp in  X_U_T   ]  

        #return super().compute_persistence(verbose, show_diag, compute_basis)
        
        if verbose:
            print('simplices',ordered_simplices)
        
        columns = [ (sp.dim, sp.boundary )    for X_U_T in self.X_U for sp in  X_U_T]  

        if verbose:
            print('boundary columns')
            for i,col in enumerate(columns):
                print(i,":",col)
        # for coloring localized cycles
        cover = [ '-'.join([str(s) for s in sp.delta ])  for X_U_T in self.X_U  for sp in  X_U_T   ]  
        if verbose:
            print('cover',cover)
        self.dgms,self.pairs  = persis.compute_persistence_dgm(ordered_simplices, 
                                               columns,cover, 
                                               show_diag=show_diag,verbose=verbose) 
        
        