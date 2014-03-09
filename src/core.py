import math
import sys
import os.path
from loader import *
from scaffold_graph import *

class Core:
	def __init__(self, clusterFile, contigFile, scaffoldFile):
		loader = Loader(clusterFile, contigFile, scaffoldFile)
		self.edgeData = loader.loadClusterFile()
		self.contigData = loader.loadContigFile()
		
		self.precomputedSums = {} #precomputed sums of contigs in the edges 
		
		self.scaffoldData, self.orderedContigsList, self.origin = loader.loadScaffoldFile()
		
		self.filterUnnecessaryEdges(5) #remove edges with less then 5 reads that supports them
		self.filterInitialUnhappyEdges() #remove edges where contigs have different orientation
		
		self.precomputeContigSums() #precompute sum of contigs inside edges
		#print self.edgeData
		
		
	def doLocalSearch(self):
		pass
		
		
	def swapContigs(self, contigID1, contigID2):
		pass
		
		
	def updateData(self):
		pass
		
		
	def filterInitialUnhappyEdges(self): #filter edges where boundary contigs have different orientation
		validOrientedEdges = []
		for edge in self.edgeData:
			if edge[1] == edge[3]:
				validOrientedEdges.append(edge)
		self.edgeData = validOrientedEdges[:]
		
		
	def getContisInBetween(self, contigID1, contigID2): #get all contigs inside edge
		firstIndex = self.orderedContigsList.index(contigID1)
		secondIndex = self.orderedContigsList.index(contigID2)
		if firstIndex > secondIndex:
			(firstIndex, secondIndex) = secondIndex, firstIndex
		contigsBetween = self.orderedContigsList[firstIndex:(secondIndex + 1)]
		return contigsBetween
		
		
	def precomputeContigSums(self): #get sum of contigs inside edge
		for edge in self.edgeData:
			contigID1 = edge[0]
			contigID2 = edge[2]
			contigsBetween = self.getContisInBetween(edge[0], edge[2])
			print contigsBetween
			currentSum = 0
			for contigID in contigsBetween:
				currentSum += self.contigData[contigID][0]
			print "sum", currentSum
			self.precomputedSums[(contigID1, contigID2)] = currentSum
		print self.precomputedSums
		return self.precomputedSums
		
		
	def precomputeGaps(self):
		pass
		
		
	def filterUnnecessaryEdges(self, size): #filter all edges with the size less then parameter size (5 in current impl.)
		validEdges = []
		for edge in self.edgeData:
			if (int(edge[-1]) > size): #check if size parameter is less then 5
				validEdges.append(edge)
		self.edgeData = validEdges[:]
		
		
	def isEdgeHappy(self, edge):
		pass
		
		
	def areContigsInSameScaffold(self, contigID1, contigID2): #we only swap contigs inside same scaffold
		if self.origin[contigID1] == self.origin[contigID2]:
			return True
		return False
		
		
		
if __name__ == "__main__":
	core = Core("drosophila_cluster.dat","drosophila.contig","drosophila.scaf")
	
