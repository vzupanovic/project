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
		
		for edge in self.edgeData: #first create empty adjency list for every node (contig)
			self.scaffoldGraph[edge[0]] = ([],[]) #first list - in edges, second list - out egdes
		
		for edge in self.edgeData: #collect all the left nodes
			wantedNode = edge[0]
			if edge[1] == '-': #in edge
				self.scaffoldGraph[wantedNode][0].append(edge[2])
			else: #out edge
				self.scaffoldGraph[wantedNode][1].append(edge[2])
			
		for edge in self.edgeData: #insert all the right nodes
			wantedNode = edge[2]
			if wantedNode in self.scaffoldGraph:
				if edge[3] == '+': #in edge 
					self.scaffoldGraph[wantedNode][0].append(edge[0])
				else: #out edge
					self.scaffoldGraph[wantedNode][1].append(edge[0])
		
		
			
			
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
		
		
	def getEdges(self, contig): #return edge data by key contig
		return self.scaffoldGraph[contig]
		
	
	def getAllPaths(self, start, end, direction = 1, path = []): #initial direction is 1 -> look out edges
		path = path + [start]
		if start == end:
			return [path]
		if not self.scaffoldGraph.has_key(start):
			return []
		paths = []
		for node in self.scaffoldGraph[start][direction]:
			direction = direction ^ 1
			if node not in path:
				newPaths = self.getAllPaths(node, end, direction, path)
				if not newPaths:
					pass
				else:
					for newPath in newPaths:
						paths.append(newPath)
		return paths
		
		
	def getShortestPath(self, start, end): #get smallest path in between two contigs
		allPaths = self.getAllPaths(start, end)
		if not allPaths:
			return None
		smallestPath = allPaths[0]
		for path in allPaths:
			if len(path) < len(smallestPath):
				smallestPath = path
		return smallestPath
		
	
	def getNextNodes(self, node): #get next node
		return self.scaffoldGraph[node][1]
			


if __name__ == "__main__":
	loader = Loader("drosophila_cluster.dat","drosophila.contig","drosophila.scaf")
	edgeData = loader.loadClusterFile()
	contigData = loader.loadContigFile()
	scaffoldData = loader.loadScaffoldFile()
	loader.loadScaffoldFile()
	scaffoldGraph = ScaffoldGraph(contigData, scaffoldData, edgeData)
		
	
	
