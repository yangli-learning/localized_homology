####################################################################################
#
#               HyperCube Cover
#
####################################################################################
#Input: 点坐标和边(with length),
#Output: 生成一组cover，X0,...,Xn
#    cover格式契合BlowupComplex的输入格式，[[point1], [point2], ..., [pointi, pointj],...]，可直接作为BlowupComplex的输入
#
####################################################################################
import random
import numpy as np
import math
import itertools
import matplotlib.pyplot as plt
import matplotlib.collections as mcoll
from pdb import set_trace as bp



class HyperCube:

    def __init__(self, point, edge, N, M, r) -> None:
      self.point = point #一组带坐标的点,numpy array
      self.edge = edge #edge_index,调用端点坐标时使用 self.point[self.edge[i][0]], self.point[self.edge[i][1]]
      self.N = N #x方向cover个数
      self.M = M #y方向cover个数
      self.r = r #0-2,overlapping degree
      self.edge_matrix= zeros_array = np.zeros((len(self.point),len(self.point)))
      self.bounder_list_x = []
      self.bounder_list_y = []
      for edge in self.edge:
        self.edge_matrix[edge[0],edge[1]] = 1
        self.edge_matrix[edge[1],edge[0]] = 1
      self.compute_cover()

    def distance_btp(p1, p2):
      r0 = pow((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2, 0.5)
      return r0

    def compute_cover(self):
      self.cover = []

      max_x = max(self.point[:,0])
      min_x = min(self.point[:,0])
      max_y = max(self.point[:,1])
      min_y = min(self.point[:,1])
      range_x = (max_x - min_x)/self.N
      range_y = (max_y - min_y)/self.M

      cover_x = []
      for i in range(self.N):
        cover_x.append([])
        self.bounder_list_x.append([min_x + i * range_x - range_x*self.r/2,
                                    min_x + (i+1)*range_x + range_x*self.r/2])
        
        for j, point in enumerate(self.point):

          if ((min_x + (i+1)*range_x + range_x*self.r/2) >= point[0] >= (min_x + i*range_x-range_x*self.r/2)):
            cover_x[i].append([j])

      cover_y = []
      for i in range(self.M):
        cover_y.append([])
        self.bounder_list_y.append([min_y + i * range_y - range_y*self.r/2,
                                    min_y + (i+1)*range_y + range_y*self.r/2])
        
        for j, point in enumerate(self.point):
          if ((min_y + (i+1)*range_y + range_y*self.r/2) >= point[1] >= (min_y + i * range_y - range_y*self.r/2)):
            cover_y[i].append([j])

      self.cover = [[] for _ in range(self.N * self.M)]
      i = 0


      for x in range(self.N):
        for y in range(self.M):
          for point in cover_x[x]:
            if point in cover_y[y]:
              self.cover[i].append(point)

          i = i+1


      step_cover = self.cover
      for j in range(len(step_cover)):

        step_edge_matrix = self.edge_matrix

        for k in range(len(step_cover[j])):

          for i in range(len(self.point)):
            
            if step_edge_matrix[step_cover[j][k][0],i] == 1 and step_edge_matrix[i,step_cover[j][k][0]] == 1 and [i] in step_cover[j]:
              if step_cover[j][k][0] < i:
                if [step_cover[j][k][0],i] in self.cover[j]:
                  continue
                else:
                  self.cover[j].append([step_cover[j][k][0],i])
              else:
                if [i, step_cover[j][k][0]] in self.cover[j]:
                  continue
                else:
                  self.cover[j].append([i,step_cover[j][k][0]])
      #return cover

    def draw_cover(self):
        for i in range(self.N):
          for j in range(self.M):
            x1, y1 = self.bounder_list_x[i][0], self.bounder_list_y[j][0]
            x2, y2 = self.bounder_list_x[i][1], self.bounder_list_y[j][1]
            plt.plot([x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], 'o-',alpha = 0.5,
                     color = [ round(i/self.N, 1),round(j/self.M, 1),round(i*j/(self.N*self.M),1)])
        
        for edge in self.edge:

            plt.plot(self.point[edge, 0], self.point[edge, 1],color= 'blue')
            #plt.triplot(self.data[:, 0], self.data[:, 1], self.simplices)
        plt.scatter(self.point[:, 0], self.point[:, 1], color='r')
        plt.show()

        return 0


    def visualize_cover_subsets(self ):
        # show each cover subset in individual subplots
        M = self.M
        N = self.N
        fig, axs = plt.subplots(M, N, figsize=(N * 4, M * 3))

        for i in range(M): #row
            for j in range(N): #col
                if M > 1 and N>1:
                    ax = axs[ i,j]
                elif N>1:
                    ax = axs[j] 
                elif M>1:
                    ax = axs[i]
                self.draw_cover_subset_by_id(i,j,ax)
        plt.tight_layout()
        plt.show()

    def draw_cover_subset_by_id(self,j,i,ax):
        x1, y1 = self.bounder_list_x[i][0], self.bounder_list_y[j][0]
        x2, y2 = self.bounder_list_x[i][1], self.bounder_list_y[j][1]
        ax.plot([x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1],
                  'ro-',alpha = 0.5)
        cover_ij = self.cover[i*self.M+j] 
        for  edge in self.edge:
          edge.sort()
          if edge in  cover_ij:
            
            ax.plot(self.point[edge, 0], self.point[edge, 1],color= 'red')
          else:
            ax.plot(self.point[edge, 0], self.point[edge, 1],color= 'gray',alpha=0.3)

        points = [ p[0] for p in cover_ij if len(p)==1]
        ax.scatter(self.point[points, 0], self.point[points, 1], color='r')

        ax.set_title("ID=%d" % (  i*self.M+j))

class HyperCubeTri(HyperCube):
  def __init__(self,point,edge,N,M,r)-> None:
      super().__init__(point,edge,N,M,r)
      self.cover =[sorted(cover, key=lambda x: (len(x ) , *x )) for cover in self.cover]


  def compute_cover(self):
    
      #self.cover = []

      max_x = max(self.point[:,0])
      min_x = min(self.point[:,0])
      max_y = max(self.point[:,1])
      min_y = min(self.point[:,1])
      range_x = (max_x - min_x)/self.N
      range_y = (max_y - min_y)/self.M

      cover_x = []
      #print(type(self.N))
      for i in range(self.N):
        cover_x.append([])
        self.bounder_list_x.append([min_x + i * range_x - range_x*self.r/2,
                                    min_x + (i+1)*range_x + range_x*self.r/2])
        #print(min_x + (i+1)*range_x + range_x*self.r , min_x + i * range_x)
        for j, point in enumerate(self.point):

          if ((min_x + (i+1)*range_x + range_x*self.r/2) >= point[0] >= (min_x + i*range_x-range_x*self.r/2)):
            cover_x[i].append([j])

      cover_y = []
      for i in range(self.M):
        cover_y.append([])
        self.bounder_list_y.append([min_y + i * range_y - range_y*self.r/2,
                                    min_y + (i+1)*range_y + range_y*self.r/2])
        #print(min_y + (i+1)*range_y + range_y*self.r , min_y + i * range_y)
        for j, point in enumerate(self.point):
          if ((min_y + (i+1)*range_y + range_y*self.r/2) >= point[1] >= (min_y + i * range_y - range_y*self.r/2)):
            cover_y[i].append([j])

      self.cover = [[] for _ in range(self.N * self.M)]
      i = 0


      for x in range(self.N):
        for y in range(self.M):
          for point in cover_x[x]:
            if point in cover_y[y]:
              self.cover[i].append(point)

          i = i+1


      step_cover = self.cover
      for j in range(len(step_cover)):

        step_edge_matrix = self.edge_matrix

        for k in range(len(step_cover[j])):

          for i in range(len(self.point)):

            #print(i,j)
            if step_edge_matrix[step_cover[j][k][0],i] == 1 and step_edge_matrix[i,step_cover[j][k][0]] == 1 and [i] in step_cover[j]:
              if step_cover[j][k][0] < i:
                if [step_cover[j][k][0],i] in self.cover[j]:
                  continue
                else:
                  self.cover[j].append([step_cover[j][k][0],i])
              else:
                if [i, step_cover[j][k][0]] in self.cover[j]:
                  continue
                else:
                  self.cover[j].append([i,step_cover[j][k][0]])

      def find_triangles(adjacency_matrix):
        num_vertices = len(adjacency_matrix)
        triangles = []
        for i in range(num_vertices):
          for j in range(i + 1, num_vertices):
            for k in range(j + 1, num_vertices):
                if adjacency_matrix[i, j] and adjacency_matrix[j, k] and adjacency_matrix[k, i]:
                    triangle = [i, j, k]
                    triangles.append(triangle)
        return triangles

      triangles_list = find_triangles(self.edge_matrix)
      step_cover = [[] for _ in range(self.N * self.M)]
      for triangle in triangles_list:
        j = -1
        for small_cover in self.cover:
          j = j+1

          if [triangle[0],triangle[1]] in small_cover and [triangle[0],triangle[2]] in small_cover and [triangle[1],triangle[2]] in small_cover:
            small_cover.append(triangle)
            continue

          if [triangle[0],triangle[1]] in small_cover or [triangle[0],triangle[2]] in small_cover or [triangle[1],triangle[2]] in small_cover:# or [triangle[0]] in small_cover  or [triangle[1]] in small_cover or [triangle[2]] in small_cover:
            if [triangle[0],triangle[1]] not in small_cover and [triangle[0],triangle[1]] not in step_cover[j]:
              step_cover[j].append([triangle[0],triangle[1]])
            if [triangle[0],triangle[2]] not in small_cover and [triangle[0],triangle[2]] not in step_cover[j]:
              step_cover[j].append([triangle[0],triangle[2]])
            if [triangle[1],triangle[2]] not in small_cover and [triangle[1],triangle[2]] not in step_cover[j]:
              step_cover[j].append([triangle[1],triangle[2]])
            if [triangle[0]] not in small_cover and [triangle[0]] not in step_cover[j]:
              step_cover[j].append([triangle[0]])
            if [triangle[1]] not in small_cover and [triangle[1]] not in step_cover[j]:
              step_cover[j].append([triangle[1]])
            if [triangle[2]] not in small_cover and [triangle[2]] not in step_cover[j]:
              step_cover[j].append([triangle[2]])
            step_cover[j].append(triangle)


      for j in range(len(step_cover)):
        for i in range(len(step_cover[j])):
          self.cover[j].append(step_cover[j][i])


  def draw_cover(self):
      #super().draw_cover(self)
      
      for i in range(self.N):
        for j in range(self.M):
          x1, y1 = self.bounder_list_x[i][0], self.bounder_list_y[j][0]
          x2, y2 = self.bounder_list_x[i][1], self.bounder_list_y[j][1]
          color_ij = [ round(i/self.N, 1),round(j/self.M, 1),round(i*j/(self.N*self.M),1)]
          plt.plot([x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], 'o-',alpha = 0.5,
                    color = color_ij)
          
          cover_ij = self.cover[i*self.M+j]  
          
          triangles = [c for c in cover_ij if len(c)==3]
          triangle_coords = [[self.point[vertex] for vertex in triangle] \
                             for triangle in triangles]

          poly_collection = mcoll.PolyCollection(triangle_coords, facecolor=color_ij,alpha=0.3)
          plt.gca().add_collection(poly_collection) 
          

      for edge in self.edge:

          plt.plot(self.point[edge, 0], self.point[edge, 1],color= 'blue')
          #plt.triplot(self.data[:, 0], self.data[:, 1], self.simplices)
      plt.scatter(self.point[:, 0], self.point[:, 1], color='r')
      plt.show()

      return 0
  
  def draw_cover_subset(self,j,i,ax):
      super().draw_cover_subset(j,i,ax)
      
      # add triangle faces in the cover 
      cover_ij = self.cover[i*self.M+j]  
      triangles = [c for c in cover_ij if len(c)==3]
      triangle_coords = [[self.point[vertex] for vertex in triangle] for triangle in triangles]

      poly_collection = mcoll.PolyCollection(triangle_coords, facecolor='red',alpha=0.5)
      ax.add_collection(poly_collection)  

      