import numpy as np
import math
import sys
import os.path
from loader import *

class ContigGraph:
	def __init__(self, orderedContigsList): #simple graph structure for contigs
	
		self.simpleContigGraph = {}
		
		sys.setrecursionlimit(100000) #sample recursion limit, have to change this
		
		copyContigList = orderedContigsList[1:]
		
		pointer = 0
		
		for contig in orderedContigsList: #get contigs from list and put them in the graph like structure
			if pointer < len(copyContigList):
				self.simpleContigGraph[contig] = copyContigList[pointer]
			else:
				self.simpleContigGraph[contig] = 'EOF'
				
			pointer += 1
				
		print self.simpleContigGraph
		
		
	def getSuccessor(self, contigID):
		return self.simpleContigGraph[contigID]
		
		
	def getPath(self, start, end, path = []): #return list of contigs between two border contigs
		path = path + [start]
		if start == end:
			return path
		if not self.simpleContigGraph.has_key(start):
			return None
		node = self.simpleContigGraph[start]
		print "tu", node
		if node not in path:
			newPath = self.getPath(node, end, path)
			if newPath:
				return newPath
		return None
		
	def swapContigs(self, contig1, contig2):
		pass
		
