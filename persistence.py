
"""
Compute the persistence homology of blowup chain complex
""" 

import numpy as np
import matplotlib.pyplot as plt
import cechmate as cm  
import phat
import itertools

def _process_distances(pairs, ordered_simplices,show_diag=False,cover=None):
    """ Setup persistence diagrams by reading off distances
    new options:
    - show_diag: keep persistence pairs with birth=death
    - cover: a list of length =len(ordered_simplices) to get the subset label of each simplex
    """ 
    dgms = {}
    # identify each simplex as a creator (1) or destroyer (-1)
    posneg = np.zeros(len(ordered_simplices))

    for [bi, di] in pairs:
        bidxs, bd = ordered_simplices[bi]  # get birth time (timestep of the creator)
        didxs, dd = ordered_simplices[di]  # get death time (timestep of the corresponding destroyer)
        assert posneg[bi] == 0 and posneg[di] == 0
        posneg[bi], posneg[di] = 1, -1     # set creator and destroyer

        assert dd >= bd
        # assert len(bidxs) == len(didxs) - 1

        p = len(bidxs) - 1    # get the dimension of this cycle (todo: use dim(sigma)+dim(Deltq) )

        if cover:
            # Don't add zero persistence pairs
            if bd != dd or  show_diag:
                dgms.setdefault(p,[]).append( [bd,dd,cover[bi]]) 
        else:
            # Don't add zero persistence pairs
            if bd != dd or  show_diag:
                dgms.setdefault(p, []).append([bd, dd])  # remove persistence pair if death = birth (on the diagonal)  
    return dgms

def _add_unpaired(dgms, pairs, simplices,cover=None):
    posneg = np.zeros(len(simplices))
    for [bi, di] in pairs:
        assert posneg[bi] == 0
        assert posneg[di] == 0
        posneg[bi] = 1
        posneg[di] = -1

    for i in range(len(posneg)):
        if posneg[i] == 0:
            (idxs, dist) = simplices[i]
            p = len(idxs) - 1
            if not p in dgms:
                dgms[p] = []
            if cover:
                dgms[p].append([dist, np.inf,cover[i]])
            else:
                dgms[p].append([dist, np.inf])

    return dgms


def max_non_infinite(numbers):
    # Filter out infinite values
    finite_numbers = [num for num in numbers if num != float('inf')  ]
    if not finite_numbers:
        return 10
    return max(finite_numbers)
    
def plotbarcode_BC(dgms,dim,coverlabel2id):
    persistence = [ [coverlabel2id[c],(b,d)]  for b,d,c in dgms[dim] ]
    labels = [ c for b,d,c in dgms[dim]]
    
    _, axes = plt.subplots(1, 1,figsize=(4,1.5))
    colormap = plt.cm.tab20.colors
    print(persistence)
    
    death = [d for b,d,c in dgms[dim]]
    
    largest_non_infinite = max_non_infinite( death )
    
    infinity  = largest_non_infinite*1.5
    
    x = [birth for (dim, (birth, death)) in persistence]
    y = [(death - birth) if death != float("inf") else (infinity - birth) for (dim, (birth, death)) in persistence]
    c = [colormap[dim % len(colormap) ] for (dim, (birth, death)) in persistence]
    
    axes.barh(range(len(x)), y, left=x, alpha=0.5, color=c, linewidth=0,label=labels)
    plt.axvline(x=infinity , color='red', linestyle='--') 


    # Shrink current axis by 20%
    box = axes.get_position()
    axes.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    axes.legend(loc='center left', bbox_to_anchor=(1, 0.5),ncol = 3)#oc="lower right")

    #axes.legend( )
    axes.set_title("Persistence barcode", fontsize=12)
    axes.set_yticks([])   #axes.invert_yaxis()
    plt.show()

def compute_persistence(ordered_simplices, columns,cover ,show_diag=False):
    if show_diag:
        print('display homology classes on the diagonal')
    

    ## Setup boundary matrix and reduce
    boundary_matrix = phat.boundary_matrix(
        columns=columns , representation=phat.representations.sparse_pivot_column
    )
    pairs = boundary_matrix.compute_persistence_pairs()
    pairs.sort()
    
    # compute persistence pairs (hide diagonals to make it easier to read)
    dgms =  _process_distances(pairs, ordered_simplices,show_diag,cover)
    for p in pairs:
        print(p)
    print(dgms)
    dgms =  _add_unpaired(dgms, pairs, ordered_simplices,cover)
    
    # assign a color id to each subset combination 
    coverlabel2id = dict([ (c,i) for i,c in enumerate(set(cover)) ])
    for dim,dgm in dgms.items(): 
        print('H%d'%dim,dgm )
        plotbarcode_BC(dgms,dim,coverlabel2id )  
    return dgms
    
