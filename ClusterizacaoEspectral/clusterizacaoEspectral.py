#!C:\python\python.exe
# coding=utf-8
from __future__ import division
from datetime import datetime
import math
import random
import argparse
import numpy
import sys
"""
This is a pure Python implementation of the K-means Clustering algorithmn. The
original can be found here:
http://pandoricweb.tumblr.com/post/8646701677/python-implementation-of-the-k-means-clustering
"""

plotly = False
try:
    import plotly
    from plotly.graph_objs import Scatter, Scatter3d, Layout
except ImportError:
    print ("INFO: Plotly is not installed, plots will not be generated.")

def main():

    ### Main program

    # Parse command line arguments in order to set simulation parameters
    # Files from simulator
    # 	groupFile.txt
    # 	coordinate.txt
    # 	distance_and_signalLoss.txt
    # 	reception.txt
   
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-n", "--numberOfSTAs", help="number of STAs in the simulation", type=int, default=1)
    parser.add_argument("-g", "--numberOfGroups", help="number of RAW grous in the simulation", type=int, default=1)
    parser.add_argument("-W", "--scenarioWidth", help="width of the area used for positioning nodes in m", type=int, default=1000)
    parser.add_argument("-H", "--scenarioHeight", help="height of the area used for positioning nodes in m", type=int, default=1000)
    parser.add_argument("-p", "--power", type=str, help="power", default=None)
    parser.add_argument("-c", "--coord", type=str, help="coord", default=None)
    global args
    args = parser.parse_args()
    
    # How many points are in our dataset?
    num_points = int(args.numberOfSTAs)

    # For each of those points how many dimensions do they have?
    # Note: Plotting will only work in two or three dimensions
    dimensions = 2   #largura e altura

    # Bounds for the values of those points in each dimension
    lower = args.scenarioWidth
    upper = args.scenarioHeight # metros

    # The K in k-means. How many clusters do we assume exist?
    #   - Must be less than num_points
    num_clusters = args.numberOfGroups

    # When do we say the process has 'converged' and stop updating clusters?
    cutoff = 0.000000000000000000000001   # 1^-25

    global eigenvectors, p_enlace, power, powerX, sortAutovals, nodes_autovec, limK, points_coord
    p_enlace = -70
    power=[]
    powerX=[]
    eigenvectors = []
    columns = []
    nodes_autovec = []
    limK = num_clusters
    points_coord = []

    # Part 1 - Spectral Algorithm
    points = SpectralAlgorithm(num_clusters, num_points, cutoff)
    
    # Part 2 - Clustering graph - K-means
    iteration_count = 20
    clusters_unbalanced = iterative_kmeans(
        points,
        num_clusters,
        cutoff,
        iteration_count
    )
    # coordinates, only for visualization
    f = open(args.coord, 'r')
    for line in f: 
        id, largura, altura = line.split(' ')
        altura = altura.replace('\n','')
        points_coord.append(Point([float(largura), float(altura)], id))
    points_coord = points_coord[0:]
    n=1
    for iz, cz in enumerate(clusters_unbalanced): 
        for pz in cz.points:
            pz.coords = points_coord[pz.coordOriginal]  #points_coord[cluster_node[n][1]]
            n = n + 1
    dimensions = 2
    # Display clusters using plotly for 2d data
    if dimensions in [2, 3] and plotly:
        plotClusters(clusters_unbalanced, dimensions)
        
    # Part 3 - Post-Processing
    list_hiddens(clusters_unbalanced) 
    clusters_balanced = post_processing(num_points, num_clusters, clusters_unbalanced)

    # generate group file for testing simulation
	# Opens GroupFile output file
    gf = open('groupFile.txt', 'w')
    for iz, cz in enumerate(clusters_balanced): 
        for pz in cz.points:
            gf.write(str(iz) + ' ' + str(pz) + '\n' )	

    gf = open('groupFilekmeans.txt', 'w')
    for iz, cz in enumerate(clusters_balanced): 
        for pz in cz.points:
            gf.write(str(iz) + ' ' + str(pz.coordOriginal) + '\n' )	

    # Display clusters using plotly for 2d data
    if dimensions in [2, 3] and plotly:
        plotClusters(clusters_balanced, dimensions)
    

def SpectralAlgorithm(num_clusters, num_points, cutoff):    
    #Generating Adjacency Matrix
    ocultoX = 0 
    m_adjacencia = numpy.zeros((num_points, num_points))
    fr = open('power.txt', 'r')
    for line in fr: 
        linha=[]
        id, a, b, c = line.split(' ')
        linha.append(Power(int(a)-1,int(b)-1,c))
        if float(c) <= float(p_enlace):
           power.append(linha)
           if float(b) > float(a):
              powerX.append(id)                      
              ocultoX = ocultoX + 1  
        m_adjacencia.itemset((int(a)-1,int(b)-1), (float(c)))  

	#Generating Degree Matrix
    m_grau = numpy.zeros((num_points, num_points))
    for x in range(num_points):
        i = 0
        for y in range(num_points):
            i = i + m_adjacencia.item((x,y)) 
        m_grau.itemset((x,x),i)    
	
    #Generating Laplacian matrix
    m_laplaciana = numpy.subtract(m_grau,m_adjacencia)
	
	#Calculating eigenvector and eigenvalue
    global autovals, autovecs, iMenor
    autovals, autovecs = numpy.linalg.eigh(m_laplaciana)

	#Eigenvalues in ascending order
    sortAutovals = autovals.copy() #numpy.array([])
    for i in range(num_points-1, 0, -1):
        for j in range(0,i):
            if sortAutovals[j] > sortAutovals[j+1]:
                swap(sortAutovals, j, j+1)
    iMenor = []
    #limK = num_clusters 
    i = 0

    while i < limK + 1:
        for j in range(num_points):
            if sortAutovals[i] == autovals[j]:
               if i > 0: # i=0 ignore first eigenvalue
				  # Select eigenvalue key
                  iMenor.append(j)
               i=i+1	
            if i > limK:
               break			   

    points = []
    commandLine = []
    posicaoOriginal = 1
    for w in range(num_points):
        # eigenvector ajusts
        var = str(autovecs[w,:]).replace('[','')  # invert
        var = str(var).replace(']','')
        var = str(var).replace('     ',',')
        var = str(var).replace('    ',',')
        var = str(var).replace('   ',',')
        var = str(var).replace('  ',',')
        var = str(var).replace(' ',',')
        if var[0] == ',':
           var = var[1:] 
        if var[len(var)-1] == ',':
           var = var.rstrip(',')
        var = str(var).replace('\n','')
        var = str(var).replace(',,',',')
        z = float(w)
        commandLine =  'points.append(Point([' + var + '], ' + str(posicaoOriginal) + '))'
        eval(commandLine)
        commandLineNodes = str(z) + ',' + var
        nodes_autovec.append(autovecs[w,1]) # invert
        posicaoOriginal = posicaoOriginal + 1
    print(points)
    

            
    return points



def iterative_kmeans(points, num_clusters, cutoff, iteration_count):
    """
    K-means isn't guaranteed to get the best answer the first time. It might
    get stuck in a "local minimum."

    Here we run kmeans() *iteration_count* times to increase the chance of
    getting a good answer.

    Returns the best set of clusters found.
    """
    candidate_clusters = []
    errors = []
    for _ in range(iteration_count):
        clusters = kmeans(points, num_clusters, cutoff)
        error = calculateError(clusters)
        candidate_clusters.append(clusters)
        errors.append(error)

    highest_error = max(errors)
    lowest_error = min(errors)
    ind_of_lowest_error = errors.index(lowest_error)
    clusters_unbalanced = candidate_clusters[ind_of_lowest_error]

    return clusters_unbalanced

def kmeans(points, k, cutoff):

    # Pick out k random points to use as our initial centroids
    initial_centroids = random.sample(points, k)

    # Create k clusters using those centroids
    # Note: Cluster takes lists, so we wrap each point in a list here.
    clusters = [Cluster([p]) for p in initial_centroids]
    # Loop through the dataset until the clusters stabilize
    loopCounter = 0
    while True:
        # Create a list of lists to hold the points in each cluster
        lists = [[] for _ in clusters]
        clusterCount = len(clusters)
        # Start counting loops
        loopCounter += 1
        # For every point in the dataset ...
        
        position_vector = 0
        for p in points:
            # Get the distance between that point and the centroid of the first
            # cluster.
            smallest_distance = getDistance(p, clusters[0].centroid, position_vector)

            # Set the cluster this point belongs to
            clusterIndex = 0

            # For the remainder of the clusters ...
            for i in range(0, clusterCount):
                # calculate the distance of that point to each other cluster's
                # centroid.
                distance = getDistance(p, clusters[i].centroid, position_vector)
                # If it's closer to that cluster's centroid update what we
                # think the smallest distance is
                if distance < smallest_distance:
                    smallest_distance = distance
                    clusterIndex = i
				   
            # After finding the cluster the smallest distance away
            # set the point to belong to that cluster
            lists[clusterIndex].append(p)
            position_vector = position_vector + 1
        # Set our biggest_shift to zero for this iteration
        biggest_shift = 0.0
        position_vector = 0

        # For each cluster ...
        for i in range(clusterCount):
            # Calculate how far the centroid moved in this iteration
            shift = clusters[i].update(lists[i])
            # Keep track of the largest move from all cluster centroid updates
            biggest_shift = max(biggest_shift, shift)

        # Remove empty clusters
        clusters = [c for c in clusters if len(c.points) != 0]
		
        # If the centroids have stopped moving much, say we're done!
        if abs(biggest_shift) < cutoff:
            break
       
    return clusters


#############################################################################
# Classes

class Point(object):
    '''
    A point in n dimensional space
    '''
    def __init__(self, coords, coordOriginal):
        '''
        coords - A list of values, one per dimension
        '''
        self.coords = coords
        self.n = len(coords)
        self.coordOriginal = coordOriginal # contains the node number - original position in the coordinate list
    def __repr__(self):
        return str(self.coords)

    def getCoord(self, n):
        return self.coords[n]

    def getCoords(self):
        return self.coords

    def getCoordOriginal(self):
        return self.coordOriginal

    def setCoord(self, n, c):
        self.coords[n] = c
        return

class Power(object):

    def __init__(self, a, b, c):
        self.a=a
        self.b=b
        self.c=c

    def getA(self):
        return self.a

    def getB(self):
        return self.b

    def getC(self):
        return self.c


class Reception(object):

    def __init__(self, a, b, c):
        self.a=a
        self.b=b
        self.c=c

    def getA(self):
        return self.a

    def getB(self):
        return self.b

    def getC(self):
        return self.c

		
class Cluster(object):
    '''
    A set of points and their centroid
    '''

    def __init__(self, points):
        '''
        points - A list of point objects
        '''

        if len(points) == 0:
            raise Exception("ERROR: empty cluster")

        # The points that belong to this cluster
        self.points = points

        # The dimensionality of the points in this cluster
        self.n = points[0].n

        # Assert that all points are of the same dimensionality
        for p in points:
            if p.n != self.n:
                raise Exception("ERROR: inconsistent dimensions")

        # Set up the initial centroid (this is usually based off one point)
        self.centroid = self.calculateCentroid()

    def __repr__(self):
        '''
        String representation of this object
        '''
        return str(self.points)

    def update(self, points):
        '''
        Returns the distance between the previous centroid and the new after
        recalculating and storing the new centroid.

        Note: Initially we expect centroids to shift around a lot and then
        gradually settle down.
        '''
        old_centroid = self.centroid
        self.points = points
        # Return early if we have no points, this cluster will get
        # cleaned up (removed) in the outer loop.
        if len(self.points) == 0:
            return 0


        self.centroid = self.calculateCentroid()
        shift = getDistanceUpdate(old_centroid, self.centroid)
        return shift

    def updatePosProcessamento(self, points):
        '''
        Returns the distance between the previous centroid and the new after
        recalculating and storing the new centroid.

        Note: Initially we expect centroids to shift around a lot and then
        gradually settle down.
        '''
        self.points = points
        # Return early if we have no points, this cluster will get
        # cleaned up (removed) in the outer loop.
        if len(self.points) == 0:
            return 0

        return 0

    def calculateCentroid(self):
        '''
        Finds a virtual center point for a group of n-dimensional points
        '''
        numPoints = len(self.points)
        # Get a list of all coordinates in this cluster
        coords = [p.coords for p in self.points]
        # Reformat that so all x's are together, all y'z etc.
        unzipped = zip(*coords)
        # Calculate the mean for each dimension
        centroid_coords = [math.fsum(dList)/numPoints for dList in unzipped]

        return Point(centroid_coords, 0)

    def getTotalDistance(self):
        '''
        Return the sum of all squared Euclidean distances between each point in 
        the cluster and the cluster's centroid.
        '''
        sumOfDistances = 0.0
        position_vector = 0
        for p in self.points:
            sumOfDistances += getDistance(p, self.centroid, position_vector)
            position_vector = position_vector + 1

        return sumOfDistances

#############################################################################
# Helper Methods

def getDistance(a, b, position_vector):
    import os
    '''
    Squared Euclidean distance between two n-dimensional points.
    https://en.wikipedia.org/wiki/Euclidean_distance#n_dimensions
    Note: This can be very slow and does not scale well
    '''
    if a.n != b.n:
        raise Exception("ERROR: non comparable points")

    accumulatedDifference = 0.0
    #call a function that creates a list with only the values of the selected eigenvectors, according to the smallest eigenvalues
    nodeAwillBeUse = True
    if nodeAwillBeUse is True:
       for i in range(a.n):
           if i > limK-1: 
              break
           squareDifference = pow((autovecs[position_vector][iMenor[i]]-b.coords[iMenor[i]]), 2)
           accumulatedDifference += squareDifference
    return accumulatedDifference

def getDistanceUpdate(a, b): 
    import os
    '''
    Squared Euclidean distance between two n-dimensional points.
    https://en.wikipedia.org/wiki/Euclidean_distance#n_dimensions
    Note: This can be very slow and does not scale well
    '''
    if a.n != b.n:
        raise Exception("ERROR: non comparable points")

    accumulatedDifference = 0.0
    for i in range(a.n):
        if i > limK-1: 
           break
        squareDifference = pow((a.coords[iMenor[i]]-b.coords[iMenor[i]]), 2)
        accumulatedDifference += squareDifference
    return accumulatedDifference

def filterA(a, q):
    c = []
    n=0
    for f in points:
        if n < q:
           m=0
           for h in autovals:
               if h == sortAutovals[n]:
                  pos = m 
                  break
               m=m+1
           c.append(points[pos])
        n=n+1
    useA = False
    for h in c:
        if h.coords == a.coords:
           useA = True 
           break
    return useA

def swap(L, i, j):
    tmp = L[i]   
    L[i] = L[j]
    L[j] = tmp
	
def makeRandomPoint(n, lower, upper):
    '''
    Returns a Point object with n dimensions and values between lower and
    upper in each of those dimensions
    '''
    p = Point([random.uniform(lower, upper) for _ in range(n)])
    return p

def calculateError(clusters):
    '''
    Return the average squared distance between each point and its cluster
    centroid.

    This is also known as the "distortion cost."
    '''
    accumulatedDistances = 0
    num_points = 0
    for cluster in clusters:
        num_points += len(cluster.points)
        accumulatedDistances += cluster.getTotalDistance()

    error = accumulatedDistances / num_points
    return error


def post_processing(num_points, num_clusters, clusters):

    for x in range(num_points):
        variavel = 'b'+str(x+1)
        globals()[variavel]=[]

    fr = open(args.power, 'r')
    for line in fr: 
        linha=[]
        id, a, b, c = line.split(' ')
        linha.append(Power(a,b,c))
        if float(c) <= float(p_enlace):
           if float(b) > 0:
              variavel = 'b' + str(b)
              globals()[variavel].append(linha)

    contaCluster = 0
    for iz, cz in enumerate(clusters):
        contaCluster = contaCluster + 1
    qtd_ideal = int(num_points/contaCluster)  #num_clusters
    looping = 0
    #print ("depois do kmeans ==================================")
    contaNodes(clusters)
    clusters = e_hiddens(clusters)
    #print ("depois e_hiddens ==================================")
    contaNodes(clusters)

    # hidden list
    print ("antes")
    list_hiddens(clusters)
	# balancing
    doador_cluster, receptor_cluster = cluster_doador_receptor(clusters, qtd_ideal)
    for receptor in receptor_cluster:    
       eleito = True  
       qtd_nos_receptor = qtdNodes(clusters, receptor)  
       while (qtd_nos_receptor < qtd_ideal):
           for doador in doador_cluster:    
              qtd_nos_doador = qtdNodes(clusters, doador)       
              if (qtd_nos_doador > qtd_ideal): 
                 # se o nó for terminal oculto retorna False
                 # se não for oculto, transfere o nó do doador para o receptor
                 eleito = foiEleito(doador, receptor, clusters, qtd_ideal)
                 if (eleito == True):
                     qtd_nos_receptor = qtd_nos_receptor + 1
                     break

           if (eleito == False):
               break
    print ("depois")
    list_hiddens(clusters)
    print ("final ==================================")
    contaNodes(clusters)
    return clusters

def foiEleito(doador, receptor, clusters, qtd_ideal):

    for izd, czd in enumerate(clusters):  # nodes of the larger cluster - doador
       if izd == doador:
          for izr, czr in enumerate(clusters): # smaller cluster - receptor
              if izr == receptor:
                 for pzd in czd.points:  # nodes of the larger cluster
                     eleito = True
                     qtd = qtdNodes(clusters, izr)
                     if qtd < (qtd_ideal):  # The cluster is below the ideal amount
                        enlace = p_enlace
                        for pzr in czr.points: # smaller cluster nodes  
                            variavel = 'b' + str(pzr.coordOriginal)
                            if float(pzr.coordOriginal) > 0:
                               nodeList = globals()[variavel]
                               for r in nodeList:  #power: 
                                   if float(pzd.coordOriginal) == float(r[0].a):
                                      if float(pzr.coordOriginal) == float(r[0].b):
                                         enlace = float(r[0].c)
                                         if enlace < p_enlace:
                                            eleito = False
                                            break

                            if enlace < p_enlace:
                               eleito = False
                               break
                        if (eleito == True):
                           czd.points.remove(pzd)
                           czr.points.append(pzd)
                           clusters[izd].updatePosProcessamento(czd.points)
                           clusters[izr].updatePosProcessamento(czr.points)
                           break  # if it was elected it leaves out of the for and if it is zero then it leaves too, to get a new larger node
          if (eleito == True):
                break
           
    return eleito
         
def contaNodes(clusters):
    print ("totais ==================================")
    for iz, cz in enumerate(clusters):
       contaNo = 0
       for pz in cz.points:
           contaNo = contaNo + 1
       print (" Cluster: ", iz, "\t -> :", contaNo)


def e_hiddens(clusters):
    oculto = 0
    for iz, cz in enumerate(clusters):
       for pz in cz.points:
           for izy, czy in enumerate(clusters):
             if iz == izy: 
               for pzy in czy.points:
                   if pz.coordOriginal != pzy.coordOriginal:
                      idX = str(pz.coordOriginal) + '->' + str(pzy.coordOriginal)
                      if idX in powerX:
                         oculto = oculto + 1
                         eleito = True
                         for outro_iz, outro_cz in enumerate(clusters):
                             if outro_iz != iz:
                                eleito = True
                                for outro_pz in outro_cz.points:
                                    idX = str(outro_pz.coordOriginal) + '->' + str(pzy.coordOriginal)
                                    if idX in powerX:
                                       oculto = oculto + 1
                                       eleito = False
                                       break
                                if eleito == True:
                                   cz.points.remove(pzy)
                                   outro_cz.points.append(pzy)
                                   clusters[izy].updatePosProcessamento(czy.points)
                                   clusters[outro_iz].updatePosProcessamento(outro_cz.points)
                                   break										   

    return clusters

def list_hiddens(clusters):
	
    oculto = 0	
    temOculto = False
    print (" Contando ocultos: ----------------------------------------------------")
    for iz, cz in enumerate(clusters):
       ocultoC = 0
       for pz in cz.points:
           for izy, czy in enumerate(clusters):
             if iz == izy: 
               for pzy in czy.points:
                   if pz.coordOriginal != pzy.coordOriginal:
                      idX = str(pz.coordOriginal) + '->' + str(pzy.coordOriginal)
                      if idX in powerX:
                         temOculto = True
                         print (" Cluster: ", iz, "\t Node :", pz.coordOriginal, pz)
           if(temOculto):
              oculto = oculto + 1
              ocultoC = ocultoC + 1
              temOculto = False
    print (" fim ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

def list_clusters(clusters):
    for iz, cz in enumerate(clusters):
       for pz in cz.points:
           for izy, czy in enumerate(clusters):
             if iz == izy: 
               for pzy in czy.points:
                   if pz.coordOriginal != pzy.coordOriginal:
                      print (" Cluster: ", iz, "\t Node :", pz.coordOriginal, pz)
								
def cont_node(clusters, cluster1, cluster2):
    conta1=0
    conta2=0
    for iz, cz in enumerate(clusters):
        for pz in cz.points:
            if (cluster1 == iz):
                conta1 = conta1 + 1
            if (cluster2 == iz):
                conta2 = conta2 + 1
    return cluster1 > cluster2

def largest_minorCluster(clusters):
    qtd_nos_receptor = 99999
    qtd_nos_doador = 0
    receptor_cluster = 99999
    doador_cluster = 0
    contador = 0
    cluster_atual = 0
    for iz, cz in enumerate(clusters):
        for pz in cz.points:
            contador = contador + 1
            if (cluster_atual != iz):
               if (contador < qtd_nos_receptor):
                    qtd_nos_receptor = contador
                    receptor_cluster = cluster_atual
               if (contador > qtd_nos_doador):
                    qtd_nos_doador = contador
                    doador_cluster = cluster_atual
               cluster_atual = iz
               contador = 0
    contador = contador + 1
    if (contador < qtd_nos_receptor):
       qtd_nos_receptor = contador
       receptor_cluster = cluster_atual
    if (contador >= qtd_nos_doador):
       qtd_nos_doador = contador
       doador_cluster = cluster_atual
      
    return qtd_nos_doador, qtd_nos_receptor 

def cluster_doador_receptor(clusters, qtd_ideal):
    #print(clusters)
    contador = 0
    doador = []
    receptor = []
    for iz, cz in enumerate(clusters):
        contador = 0
        for pz in cz.points:
            contador = contador + 1
        if (contador < qtd_ideal):
            receptor.append(iz)
        if (contador > qtd_ideal):
            doador.append(iz)
 
    return doador, receptor 

def qtdNodes(clusters, cluster):
    nos_cluster = 0
    contador = 0
    for iz, cz in enumerate(clusters):
        for pz in cz.points:
            if (cluster == iz):
               contador = contador + 1
    return contador 


def plotClusters(data, dimensions):
    '''
    This uses the plotly offline mode to create a local HTML file.
    This should open your default web browser.
    '''
    if dimensions not in [2, 3]:
        raise Exception("Plots are only available for 2 and 3 dimensional data")

    #From itertools import chain
    #Convert data into plotly format.
    traceList = []
    for i, c in enumerate(data):
        # Get a list of x,y coordinates for the points in this cluster.
        cluster_data = []
        zz = []
        for z in c.points:
            #cluster_data.append(z.coords)
            zz.append(z.getCoords())
            cluster_data = zz #z.coords

        #cluster_data = list(chain(*cluster_data))
        trace = {}
        itemTuple = ()
        centroid = {}
        if dimensions == 2:
            # Convert our list of x,y's into an x list and a y list.
            ij = []
            for j in cluster_data:
                ij.append(j.coords) 
            trace['x'], trace['y'] = zip(*ij)
            trace['mode'] = 'markers'
            trace['marker'] = {}
            trace['marker']['symbol'] = i
            trace['marker']['size'] = 12
            trace['name'] = "Cluster " + str(i)
            traceList.append(Scatter(**trace))
            traceList.append(Scatter(**centroid))
        else:
            symbols = [
                "circle",
                "square",
                "diamond",
                "circle-open",
                "square-open",
                "diamond-open",
                "cross", "x"
            ]
            symbol_count = len(symbols)
            if i > symbol_count:
                print ("Warning: Not enough marker symbols to go around")
            # Convert our list of x,y,z's separate lists.
            trace['x'], trace['y'], trace['z'] = zip(*cluster_data)
            trace['mode'] = 'markers'
            trace['marker'] = {}
            trace['marker']['symbol'] = symbols[i]
            trace['marker']['size'] = 12
            trace['name'] = "Cluster " + str(i)
            traceList.append(Scatter3d(**trace))
            # Centroid (A trace of length 1)
            centroid['x'] = [c.centroid.coords[0]]
            centroid['y'] = [c.centroid.coords[1]]
            centroid['z'] = [c.centroid.coords[2]]
            centroid['mode'] = 'markers'
            centroid['marker'] = {}
            centroid['marker']['symbol'] = symbols[i]
            centroid['marker']['color'] = 'rgb(200,10,10)'
            centroid['name'] = "Centroid " + str(i)
            traceList.append(Scatter3d(**centroid))

    title = "K-means clustering with %s clusters and %s limK" % (str(len(data)), str(limK))
    plotly.offline.plot({
        "data": traceList,
        "layout": Layout(title=title)
    })

if __name__ == "__main__":
    main()