import itertools
from sklearn import datasets
from sklearn.datasets import make_circles
import numpy as np
from collections import Counter


#from matplotlib import pyplot as plt

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

def legend_without_duplicate_labels(ax):
    # https://stackoverflow.com/questions/19385639/duplicate-items-in-legend-in-matplotlib
    handles, labels = ax.get_legend_handles_labels()
    handles.reverse()
    labels.reverse()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles , labels )) if l not in labels[:i]]
    
    ax.legend(*zip(*unique),loc='center left', 
              bbox_to_anchor=(1, 0.5), ncol = 3) 

def find_superset( sets, query):
       
    # Extract tuple supersets from List
    # Using all() + list comprehension + Counter
    res = [sub for sub in sets if
        all(Counter(sub)[x] >= Counter(query)[x]
        for x in Counter(query))]
    return res 

def transform_img2pc(img):
    img_array = np.asarray(img)
    indices = np.argwhere(img_array > 120)
    return indices.astype(np.float32)
