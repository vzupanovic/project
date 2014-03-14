import numpy as np
import math
import sys
import os.path
from loader import *

class ContigGraph:
	def __init__(self, orderedContigsList, contigInfo): #simple graph structure for contigs
	
		self.simpleContigGraph = {}
		
		sys.setrecursionlimit(100000) #sample recursion limit, have to change this
		
		copyContigList = orderedContigsList[1:]
		
		pointer = 0
		
		for contig in orderedContigsList:  #initialzie empty simple graph strucure (similar to double linked list) with gap values inside dictionary
			self.simpleContigGraph[contig] = [{},{}]
		
		for contig in orderedContigsList: #get contigs from list and put them in the graph like structure
			if pointer < len(copyContigList):
				self.simpleContigGraph[contig][1][copyContigList[pointer]] = 0
				self.simpleContigGraph[contig][1][copyContigList[pointer]] =  contigInfo[self.simpleContigGraph[contig][1].keys()[0]][2]
				
			pointer += 1
			
		for pointer in range(len(orderedContigsList)-1, 0, -1):
			contig = orderedContigsList[pointer]
			self.simpleContigGraph[contig][0][orderedContigsList[pointer - 1]] = 0
			self.simpleContigGraph[contig][0][orderedContigsList[pointer - 1]] =  contigInfo[self.simpleContigGraph[contig][0].keys()[0]][2]
		
		
		
	def getNextNode(self, contigID): #get next node
		nextKeys = self.simpleContigGraph[contigID][1].keys()
		if nextKeys:
			return nextKeys[0]
		return 'EOF'
		
	
	def getPreviousNode(self, contigID): #get previous node
		previousKeys = self.simpleContigGraph[contigID][0].keys()
		if previousKeys:
			return previousKeys[0]
		return 'BEG'
		
		
	def isPreviousNode(self, contigID, previous): #check if node is previous
		if self.getPreviousNode(contigID) == previous:
			return True
			
			
	def isNextNode(self, contigID, nextNode): #check if node is next
		if self.getNextNode(contigID) == nextNode:
			return True
		
		
	def getPath(self, start, end, direction = 1, path = []): #return list of contigs between two border contigs direction (1-front), 0-backwards
		path = path + [start]
		if start == end:
			return path
		if not self.simpleContigGraph.has_key(start):
			return None
		node = self.simpleContigGraph[start][direction].keys()[0]
		if node not in path:
			newPath = self.getPath(node, end, direction, path)
			if newPath:
				return newPath
		return None
		
	
	def getSum(self, start, end, direction = 1, contigSum = 0): #return sum of contigs between the two border contigs
		if start == end:
			return contigSum + int(contigInfo[start][1])
		if not self.simpleContigGraph.has_key(start):
			return 0
		node = self.simpleContigGraph[start][direction].keys()[0]
		if node:
			contigSum = int(contigInfo[start][1]) + self.getSum(node, end, direction, contigSum)
		return contigSum

		
		
	def basicSwap(self, contig1, contig2):
			firstChange = self.simpleContigGraph[contig2][0].keys()[0]
			lastChange = self.simpleContigGraph[contig1][1].keys()[0]
			self.simpleContigGraph[firstChange][1] = self.simpleContigGraph[contig2][1]
			self.simpleContigGraph[lastChange][0] = self.simpleContigGraph[contig1][0]
			pom = self.simpleContigGraph[contig2][0]
			self.simpleContigGraph[contig2][0] = self.simpleContigGraph[contig2][1]
			self.simpleContigGraph[contig2][1] = self.simpleContigGraph[contig1][1]
			pom1 = self.simpleContigGraph[contig1][0]
			self.simpleContigGraph[contig1][0] = pom
			self.simpleContigGraph[contig1][1] = pom1
			
			
		
	def swapContigs(self, contig1, contig2): #swap to contigs, two basic cases
		if self.isPreviousNode(contig1, contig2): #contig2 previous
			self.basicSwap(contig1, contig2)
	
			
		elif self.isPreviousNode(contig2, contig1):
			(contig2, contig1) = (contig1, contig2)
			self.basicSwap(contig1, contig2)

		else:
			firstChange1 = self.simpleContigGraph[contig1][0].keys()[0] #memorize first contig to change
			lastChange1 = self.simpleContigGraph[contig1][1].keys()[0]
			
			firstChange2 = self.simpleContigGraph[contig2][0].keys()[0] #memorize last contig to change
			lastChange2 = self.simpleContigGraph[contig2][1].keys()[0]
			
			pomNext = self.simpleContigGraph[firstChange1][1] #change next nodes of previous nodes
			self.simpleContigGraph[firstChange1][1] = self.simpleContigGraph[firstChange2][1]
			self.simpleContigGraph[firstChange2][1] = pomNext
			
			pomPrevious = self.simpleContigGraph[lastChange1][0] #change previous nodes of next nodes
			self.simpleContigGraph[lastChange1][0] = self.simpleContigGraph[lastChange2][0]
			self.simpleContigGraph[lastChange2][0] = pomPrevious
			
			pomPrevious = self.simpleContigGraph[contig1][0] #swap contigs
			pomNext = self.simpleContigGraph[contig1][1]
			self.simpleContigGraph[contig1][0] = self.simpleContigGraph[contig2][0]
			self.simpleContigGraph[contig2][0] = pomPrevious
			self.simpleContigGraph[contig1][1] = self.simpleContigGraph[contig2][1]
			self.simpleContigGraph[contig2][1] = pomNext
			
			
			
		
if __name__ == "__main__":
	loader = Loader("drosophila_cluster.dat","drosophila.contig","drosophila.scaf")
	scaffoldData, orderedContigsList, origin, contigInfo = loader.loadScaffoldFile()
	contiger = ContigGraph(orderedContigsList, contigInfo)
