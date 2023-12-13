import itertools
from sklearn import datasets
from sklearn.datasets import make_circles
import numpy as np
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

def threecircles(N,s):

    circle1 = datasets.make_circles(n_samples=N,random_state=s)[0]+ 5 * datasets.make_circles(n_samples=N,random_state=s)[0]

    circle2 = datasets.make_circles(n_samples=N,random_state=s)[0]+ 5 * datasets.make_circles(n_samples=N,random_state=s)[0]
    circle3 = datasets.make_circles(n_samples=N,random_state=s)[0]+ 5 * datasets.make_circles(n_samples=N,random_state=s)[0]
    for circle in circle2:
        circle[0] = circle[0] + 10
    for circle in circle3:
        circle[0] = circle[0] + 20

    circles_data = np.concatenate([circle1, circle2, circle3])
    return circles_data
