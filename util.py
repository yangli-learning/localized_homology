import itertools
def make_standard_simplex(n):
    # standard simplex Delta^J 
    # Delta^J  = [ ([0],0), ([1],0),([0,1],1) ]
    D_J = [] 
    for i in range(n):
        c = list(itertools.combinations(range(n),i+1)) 
        D_J = D_J + [   list(s)   for s in c]  
    return D_J

def subset_intersection(list_of_complexes): 
    # compute intersection of a list of simplicial covers
    sets = []
    for CP in list_of_complexes: 
        sets.append(  set([tuple(sigma) for sigma in CP]))
    intersection = set.intersection(*sets) 
    # sort by simplex dimension, then by simplex id
    intersection  = sorted(intersection, key=lambda x: (len(x ) , x[0] )) 
    
    return [list(t) for t in intersection]
