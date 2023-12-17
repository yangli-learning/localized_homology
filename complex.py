import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

"""
Input: np array of coordinate : [[1.1, 2.2],[2.3, 3.5]] 2D for now

Output: edge coordinates and index

Example:

complex = DelaunayComplex(np.array(data),r) # r is maximum edge length
edge_list_position = complex.edge() 
complex.edge_index()  
complex.draw_complex()  
"""
class DelaunayComplex:

    def __init__(self, data, max_r) -> None:
        # lets data X a set finite in R^d
        self.data = data
        self.simplices = Delaunay(self.data).simplices
        self.max_r = max_r

    def mesh(self):
      return self.simplices



    def edge_index(self):
      def distance_btp(p1, p2):
        r0 = pow((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2, 0.5)
        return r0

      #simplices = Delaunay(self.data).simplices
      a = self.simplices
      #print(a)
      n = len(self.data.tolist())

      edge_matrix = zeros_array = np.zeros((n, n))

      for triangle in self.simplices:
        edge_matrix[triangle[0], triangle[1]] = 1
        edge_matrix[triangle[0], triangle[2]] = 1
        edge_matrix[triangle[1], triangle[2]] = 1
        edge_matrix[triangle[1], triangle[0]] = 1
        edge_matrix[triangle[2], triangle[0]] = 1
        edge_matrix[triangle[2], triangle[1]] = 1

      #print(edge_matrix)
      edge_list_index = []

      for  i in range(n):
        for j in range(i):
          if edge_matrix[i, j] == 1:
            #print(self.data[i],self.data[j])
            if(distance_btp(self.data[i],self.data[j]) < self.max_r):
              edge_list_index.append([i,j])


      return edge_list_index
    def draw_complex(self):
      for edge in self.edge_index():
        plt.plot(self.data[edge, 0], self.data[edge, 1],color= 'blue')
      #plt.triplot(self.data[:, 0], self.data[:, 1], self.simplices)
      plt.scatter(self.data[:, 0], self.data[:, 1], color='r')
      plt.show()