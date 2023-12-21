
"""
Compute the persistence homology of blowup chain complex
""" 

import numpy as np
import matplotlib.pyplot as plt
import cechmate as cm  
import phat
import itertools
import util


def _get_barcode_from_pair(creator,destroyer, ordered_simplices):
    creator_simplex, birth= ordered_simplices[creator]
    destroyer_simplex, death = ordered_simplices[destroyer]
    if destroyer == -1: # for compatibility sometimes inf is represented by -1
        death = np.inf 
    dim = len(creator_simplex)-1 # dimension of the cycle
    return (birth,death) ,dim

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
        (bd,dd),p = _get_barcode_from_pair(bi,di,ordered_simplices)
        #bidxs, bd = ordered_simplices[bi]  # get birth time (timestep of the creator)
        #didxs, dd = ordered_simplices[di]  # get death time (timestep of the corresponding destroyer)
        assert posneg[bi] == 0 and posneg[di] == 0
        posneg[bi], posneg[di] = 1, -1     # set creator and destroyer

        assert dd >= bd
        # assert len(bidxs) == len(didxs) - 1

        #p = len(bidxs) - 1    # get the dimension of this cycle (todo: use dim(sigma)+dim(Deltq) )

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
            (idxs, dist) = simplices[i] # index is the simplex, # dist is the filtration parameter
            p = len(idxs) - 1
            if not p in dgms:
                dgms[p] = []
            if cover:
                dgms[p].append([dist, np.inf,cover[i]])
            else:
                dgms[p].append([dist, np.inf])
            pairs.append_pair(  i,  -1 )  # add infinity as "-1"
    return dgms,pairs


def _max_non_infinite(numbers):
    # Filter out infinite values
    finite_numbers = [num for num in numbers if num != float('inf')  ]
    if not finite_numbers:
        return 10 
    return max(finite_numbers)
    
def plotbarcode_BC(dgms,dim,coverlabel2id=None):
    if coverlabel2id:
        persistence = [ [coverlabel2id[c],(b,d)]  for b,d,c in dgms[dim] ]
        labels = [ c for b,d,c in dgms[dim]]
    else:
        persistence =   [ [0,(b,d)]  for b,d  in dgms[dim] ]
        labels =  ''

    print("number of bars",len(persistence))
    _, axes = plt.subplots(1, 1,figsize=(4,1.5))
    colormap = plt.cm.tab20.colors

    death = [bd[1] for bd in dgms[dim]]
    
    largest_non_infinite = _max_non_infinite( death )
    # if max non infinite is zero, 
    if largest_non_infinite == 0:
        infinity = 10
    else:
        infinity  = largest_non_infinite*1.8
    
    x = [birth for (cid, (birth, death)) in persistence]
    y = [(death - birth) if death != float("inf") else (infinity - birth) for (cid, (birth, death)) in persistence]
    if coverlabel2id:
        c = [colormap[cid % len(colormap) ] for (cid, (birth, death)) in persistence]
    else:
        c = colormap[0]
    axes.barh(range(len(x)), y, left=x, alpha=0.5, color=c, linewidth=0,label=labels)
    axes.axvline(x=infinity , color='red', linestyle='--') 


    # Shrink current axis by 20%
    box = axes.get_position()
    axes.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    util.legend_without_duplicate_labels(axes)
    
    #axes.legend( )
    axes.set_title("Persistence barcode (H%d)" % dim, fontsize=12)
    axes.set_yticks([])   #axes.invert_yaxis()
    plt.show()

def compute_persistence_dgm(ordered_simplices, columns,cover ,show_diag=False,verbose=False):
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
    dgms ,pairs =  _add_unpaired(dgms, pairs, ordered_simplices,cover)
    
    if verbose:
        for b,d in pairs:
         print(b,d )

    # assign a color id to each subset combination 
    coverlabel2id = dict([ (c,i) for i,c in enumerate(set(cover)) ])
    for dim,dgm in dgms.items(): 
        print('H%d |'%dim,dgm )
        plotbarcode_BC(dgms,dim,coverlabel2id )  
    
    return dgms, pairs  

def compute_basis_from_persistence_pairs(columns, ordered_simplices,pairs):
    cycle_basis = dict()
    for p in pairs: 
        if p[1]==-1: # only compute basis living at time n-1
            basis = find_basis(p[0],columns,pairs)
            barcode,dim = _get_barcode_from_pair(p[0],p[1], ordered_simplices)
            cycle_basis[p[0]] =  dict({ 'barcode':barcode, 'H_dim':dim,'basis':basis})
    return cycle_basis


def find_basis(sid,columns,pairs ):
    # sid is the id of a creator. we repeatedly find the cascade of sid, 
    # which forms the basis of the homology class
    cascade = set([sid])
    basis = set([sid])
    partner = [0]* len(columns) 
    for birth,death in pairs:
        partner[birth] = death
        partner[death] = birth 
    while cascade :
        sigma = cascade.pop()
        dim, boundary = columns[sigma]
        for tau in boundary:
            
            if partner[tau] != -1 and partner[tau] not in basis:
                cascade.add(partner[tau]) 
        basis.add(sigma)
    return basis


    
