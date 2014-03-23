import numpy as np
import math
import sys
import os.path
from loader import *
from scaffold import *

class ScaffoldGraph:
	def __init__(self, contigData, scaffoldData, edgeData): #contigData - dict, edgeData - list, scaffoldData - dict
		self.contigData = contigData
		self.scaffoldData = scaffoldData
		self.edgeData = edgeData
		
		self.scaffoldGraph = {} #initialize empty dict for storing graph like structure

		for edge in self.edgeData: #first create empty adjency list for every node (contig)
			self.scaffoldGraph[edge[0]] = ({},{}) #first dict - in edges, second dict - out egdes
		
		for edge in self.edgeData: #collect all the left nodes
			wantedNode = edge[0]
			if edge[1] == '-': #in edge
				self.scaffoldGraph[wantedNode][0][edge[2]] = (edge[4],edge[5],edge[6])
			else: #out edge
				self.scaffoldGraph[wantedNode][1][edge[2]] = (edge[4],edge[5],edge[6])
			
		for edge in self.edgeData: #insert all the right nodes
			wantedNode = edge[2]
			if wantedNode in self.scaffoldGraph:
				if edge[3] == '+': #in edge 
					self.scaffoldGraph[wantedNode][0][edge[0]] = (edge[4],edge[5],edge[6])
				else: #out edge
					self.scaffoldGraph[wantedNode][1][edge[0]] = (edge[4],edge[5],edge[6])
		
		
		
		
		
	
	def getDistance(self, firstContig, secondContig): #get distance between contigs (edge), this is mean!
		distancePrevious = 0
		distanceNext = 0
		if firstContig in self.scaffoldGraph:
			dataPrevious = self.scaffoldGraph[firstContig][0]
			dataNext = self.scaffoldGraph[firstContig][0]
			if secondContig in dataPrevious:
				distancePrevious = dataPrevious[secondContig][0]
			if secondContig in dataNext:
				distanceNext = dataNext[secondContig][0]
		return (distancePrevious, distanceNext)
		
		
	def getDeviation(self, firstContig, secondContig): #get standard deviation of the edge
		deviationPrevious = 0
		deviationNext = 0
		if firstContig in self.scaffoldGraph:
			dataPrevious = self.scaffoldGraph[firstContig][0]
			dataNext = self.scaffoldGraph[firstContig][0]
			if secondContig in dataPrevious:
				deviationPrevious = dataPrevious[secondContig][1]
			if secondContig in dataNext:
				deviationNext = dataNext[secondContig][1]
		return (deviationPrevious, deviationNext)
	
		
	def getSize(self, firstContig, secondContig): #get distance between two border contigs in the edge
		sizePrevious = 0
		sizeNext = 0
		if firstContig in self.scaffoldGraph:
			dataPrevious = self.scaffoldGraph[firstContig][0]
			dataNext = self.scaffoldGraph[firstContig][0]
			if secondContig in dataPrevious:
				sizePrevious = dataPrevious[secondContig][2]
			if secondContig in dataNext:
				sizeNext = dataNext[secondContig][2]
		return (sizePrevious, sizeNext)
			
		
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
		
	
	def getNextNodes(self, node): #get next nodes
		return self.scaffoldGraph[node][1]
		
		
	def getPreviousNodes(self, node): #get previous nodes
		return self.scaffoldGraph[node][0]
			


if __name__ == "__main__":
	loader = Loader("drosophila_cluster.dat","drosophila.contig","drosophila.scaf")
	edgeData, edgeLinkers = loader.loadClusterFile()
	contigData = loader.loadContigFile()
	scaffoldData = loader.loadScaffoldFile()
	loader.loadScaffoldFile()
	scaffoldGraph = ScaffoldGraph(contigData, scaffoldData, edgeData)
		
	
	
