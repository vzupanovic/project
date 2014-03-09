import numpy as np
import math
import sys
import os.path
from loader import *

class ScaffoldGraph:
	def __init__(self, contigData, scaffoldData, edgeData): #contigData - dict, edgeData - list, scaffoldData - dict
		self.contigData = contigData
		self.scaffoldData = scaffoldData
		self.edgeData = edgeData
		
		self.scaffoldGraph = {}
		
		self.filterUnnecessaryEdges(5) # edge must be supported with more then 5 reads
		
		for edge in self.edgeData:
			edgeKey = (edge[0], edge[2]) #contigs that determin the edge
			edgeDict = {} #every edge has its own dict

			edgeDict[edgeKey] = ([edge[1], edge[3], edge[4], edge[5], edge[6]],0,0) #last parameters - sum of contigs in edge, sum of gaps
			self.scaffoldGraph[edge[0]] = edgeDict
			self.scaffoldGraph[edge[1]] = edgeDict
		print self.scaffoldGraph
			
			
	def filterUnnecessaryEdges(self, size): #filter all edges with the size less then parameter size (5 in current impl.)
		validEdges = []
		for edge in self.edgeData:
			if (int(edge[-1]) > 5): #check if size parameter is less then 5
				validEdges.append(edge)
		self.edgeData = validEdges[:]
		
		
	def getContigSumInEdge(self):
		pass
		
		
	def updateSize(self, edge):
		pass
		
		
	def getEdge(self, contig): #return edge data by key contig
		return self.scaffoldGraph[contig]
		
	
			


if __name__ == "__main__":
	loader = Loader("drosophila_cluster.dat","drosophila.contig","drosophila.scaf")
	edgeData = loader.loadClusterFile()
	contigData = loader.loadContigFile()
	scaffoldData = loader.loadScaffoldFile()
	loader.loadScaffoldFile()
	scaffoldGraph = ScaffoldGraph(contigData, scaffoldData, edgeData)
		
	
	
