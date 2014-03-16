import math
import sys
import os.path
from loader import *
from scaffold_graph import *
from contig_graph import *

class Core:
	def __init__(self, clusterFile, contigFile, scaffoldFile):
		loader = Loader(clusterFile, contigFile, scaffoldFile)
		self.edgeData, self.edgeLinkers = loader.loadClusterFile()
		self.contigData = loader.loadContigFile()
		self.scaffoldData, self.orderedContigsList, self.origin, self.contigInfo = loader.loadScaffoldFile()
		
		self.precomputedSums = {} #precomputed sums of contigs in the edges 
		
		self.filterUnnecessaryEdges(5) #remove edges with less then 5 reads that supports them
		self.filterInitialUnhappyEdges() #remove edges where contigs have different orientation
		
		self.contigGraphObject = ContigGraph(self.orderedContigsList, self.contigInfo) #define simple contig graph structure
		self.scaffoldGraphObject = ScaffoldGraph(self.contigData, self.scaffoldData, self.edgeData) #define scaffold graph structure
		
		self.precomputeContigSums() #precompute sum of contigs inside edges
		#to do precompute gaps
		
		
		
		
	def doLocalSearch(self):
		pass
		
		
	def updateData(self):
		#compute gaps again??
		pass
		
		
	def filterInitialUnhappyEdges(self): #filter edges where boundary contigs have different orientation
		validOrientedEdges = []
		for edge in self.edgeData:
			if edge[1] == edge[3]:
				validOrientedEdges.append(edge)
		self.edgeData = validOrientedEdges[:]
		
		
	def precomputeContigSums(self): #get sum of contigs inside edge and add that to scaffold graph structure
		for edge in self.edgeData: #collect all the left nodes
			wantedNode = edge[0]
			direction = self.contigGraphObject.getDirection(edge[0],edge[2])
			cSum = self.contigGraphObject.getSum(edge[0], edge[2],direction)
			if edge[1] == '-': #in edge
				data = self.scaffoldGraphObject.scaffoldGraph[wantedNode][0][edge[2]]
				data = data + (cSum,)
				self.scaffoldGraphObject.scaffoldGraph[wantedNode][0][edge[2]] = data
			else: #out edge
				data = self.scaffoldGraphObject.scaffoldGraph[wantedNode][1][edge[2]]
				data = data + (cSum,)
				self.scaffoldGraphObject.scaffoldGraph[wantedNode][1][edge[2]] = data
			
		for edge in self.edgeData: #insert all the right nodes
			wantedNode = edge[2]
			if wantedNode in self.scaffoldGraphObject.scaffoldGraph:
				if edge[3] == '+': #in edge 
					data = self.scaffoldGraphObject.scaffoldGraph[wantedNode][0][edge[0]]
					data = data + (cSum,)
					self.scaffoldGraphObject.scaffoldGraph[wantedNode][0][edge[0]] = data
				else: 
					data = self.scaffoldGraphObject.scaffoldGraph[wantedNode][1][edge[0]]
					data = data + (cSum,)
					self.scaffoldGraphObject.scaffoldGraph[wantedNode][1][edge[0]] = data
		
		
			
		
		
	def precomputeGaps(self):
		pass
		
		
	def filterUnnecessaryEdges(self, size): #filter all edges with the size less then parameter size (5 in current impl.)
		validEdges = []
		for edge in self.edgeData:
			if (int(edge[-1]) > size): #check if size parameter is less then 5
				validEdges.append(edge)
		self.edgeData = validEdges[:]
		
		
	def isEdgeHappy(self, edge):
		t = 6 #must be 6, don't know why
		direction = self.contigGraphObject.getDirection(edge[0], edge[2])
		contigSum = self.contigGraphObject.getSum(edge[0], edge[2], direction)
		mean = self.scaffoldGraphObject.getDistance(edge[0], edge[2])[0]
		sigma = self.scaffoldGraphObject.getDeviation(edge[0], edge[2])[0]
		if contigSum < (mean + t*sigma):
			return True
		return False
		
		
	def areContigsInSameScaffold(self, contigID1, contigID2): #we only swap contigs inside same scaffold
		if self.origin[contigID1] == self.origin[contigID2]:
			return True
		return False
		
		
		
if __name__ == "__main__":
	core = Core("drosophila_cluster.dat","drosophila.contig","drosophila.scaf")
	
