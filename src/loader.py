###################################################################
# module for loading input files into the simple data structures
# edgeData - simple list, contigs - dictionary key is the contig name
# value is tuple (contig len, sequence), scaffoldData - dict, key
# is scaffold id (ie. opera_scaffold_319) and value is list of
# contig names, orientation, contig length and gap size.

import numpy as np
import math
import sys
import os.path
from Bio import SeqIO


class Loader: #class for loading files in simple datastructures
	def __init__(self, clusterFile, contigFile, scaffoldFile):
		self.clusterFile = clusterFile
		self.scaffoldFile = scaffoldFile
		self.contigFile = contigFile
		
		if os.path.isfile(self.clusterFile)!=True:
			print "File "+self.clusterFile+" doesn't exist, exiting..."
			exit(-1)
		if os.path.isfile(self.scaffoldFile)!=True:
			print "File "+self.scaffoldFile+" doesn't exist, exiting..."
			exit(-1)
		if os.path.isfile(self.contigFile)!=True:
			print "File "+self.contigFile+" doesn't exist, exiting..."
			exit(-1)
		
		
	def loadClusterFile(self): #load cluster file into list and return it
		self.edgeData = []
		self.clusterDataDouble = {} #dict keys are pairs of contigs, necessary for constructing scaffold graph
		stream = open(self.clusterFile, 'r')
		data = stream.readlines()
		for line in data:
			line = line.strip()
			temp = line.split('\t')
			self.edgeData.append(temp)
		self.edgeData = self.edgeData[1:]
		for edge in self.edgeData:
			self.clusterDataDouble[(edge[0],edge[2])] = edge
		return self.edgeData, self.clusterDataDouble
		
		
	def loadContigFile(self): #load contig file into dictionary, key is contig name and return it
		self.contigs = {}
		for seqRecord in SeqIO.parse(self.contigFile, "fasta"):
			self.contigs[seqRecord.id] = (len(seqRecord.seq), seqRecord.seq) #value is tuple contig len and seq. itself
		return self.contigs
	
		
	def loadScaffoldFile(self): #load scaffold file into dictionary, key is scaffold name, value - ordered list of contigs
		self.scaffoldData = {}
		self.orderedContigs = [] #list of all oredred contigs, neccessary for swaping contigs
		self.origin = {} #dict of origin, every contigs (key) has its origin scaffold
		self.contigInfo = {} #dict, key is contig name values are (orientation, len, gap size)
		stream = open(self.scaffoldFile, 'r')
		data = stream.readlines()
		for line in data:
			line = line.strip()
			if (line[0] == '>'):
				line = line[1:]
				temp = line.split("\t")
				currentHeader = temp[0]
				self.scaffoldData[currentHeader] = []
			else:
				temp = line.split('\t')
				self.scaffoldData[currentHeader].append(temp)
				self.orderedContigs.append(temp[0])
				self.origin[temp[0]] = currentHeader
				self.contigInfo[temp[0]] = (temp[1], temp[2], temp[3]) 
		return self.scaffoldData, self.orderedContigs, self.origin, self.contigInfo
		
		
