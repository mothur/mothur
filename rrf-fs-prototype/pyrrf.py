#!/usr/bin/python
from __future__ import division
from idlelib.TreeWidget import TreeNode
from math import log, ceil
import random
import inspect


# The main algorithm resides here, it the manipulator class
class DecisionTree:
	def __init__(self, baseDataSet):
		self.baseDataSet = baseDataSet
		self.numSamples = len(baseDataSet)
		self.numFeatures = len(baseDataSet[0]) - 1
		self.numOutputClasses = 0
		self.classes = []
		self.bootstrappedTrainingSamples = []
		self.bootstrappedTestSamples = []
		self.rootNode = None

		for i in range(0, self.numSamples):
			outcome = baseDataSet[i][-1]
			if outcome not in self.classes:
				self.classes.append(outcome)
				self.numOutputClasses += 1

		self.createBootStrappedSamples()
		self.buildDecisionTree()

	def createBootStrappedSamples(self):
		isInTrainingSamples = [False for i in range(0, self.numSamples)]
		for i in range(0, self.numSamples):
			index = random.randint(0, self.numSamples - 1)
			self.bootstrappedTrainingSamples.append(self.baseDataSet[index])
			isInTrainingSamples[index] = True

		for i in range(0, self.numSamples):
			if not isInTrainingSamples[i]:
				self.bootstrappedTestSamples.append(self.baseDataSet[i])
#		for x in self.bootstrappedTrainingSamples: print x
#		print "###"
#		for x in self.bootstrappedTestSamples: print x

	# randomly select log2(totalFeatures) number features for each node
	def selectFeatureSubsetRandomly(self):
		# TODO optimum feature is log2(totalFeatures), we might need to modify this one
		self.optimumFeatureSubsetSize = int(ceil(log(self.numFeatures, 2)))

		featureSubsetIndices = []
		while len(featureSubsetIndices) < self.optimumFeatureSubsetSize:
			randomIndex = random.randint(0, self.numFeatures - 1)
			if randomIndex not in featureSubsetIndices: featureSubsetIndices.append(randomIndex)
		return sorted(featureSubsetIndices)

	def buildDecisionTree(self):
		print "buildDecisionTree()"
		treeNode = TreeNode(self.bootstrappedTrainingSamples, self.numFeatures, self.numSamples, self.numOutputClasses, 0)
		self.splitRecursively(treeNode)
		self.printTree(treeNode, "root")

	def printTree(self, treeNode, caption):
		tabs = ""
		for i in range(0, treeNode.generation): tabs += '\t'

		if not treeNode.isLeaf:
			print tabs + caption + ' ( ' + str(treeNode.splitFeatureValue) + ' < X'+ str(treeNode.splitFeatureIndex) + ' )'
			self.printTree(treeNode.leftChildNode, 'leftChild')
			self.printTree(treeNode.rightChildNode, 'rightChild')
		else:
			print tabs + caption + ' ( classified to: ' + str(treeNode.outputClass) +', samples: ' + str(treeNode.numSamples) +' )'



	def splitRecursively(self, treeNode):
		print "splitRecursively()"

		# return immediately if this node is just a leaf, no recursion is needed
		if treeNode.numSamples < 2:
			print "Already classified: Case 1"
			treeNode.isLeaf = True
			treeNode.outputClass = treeNode.bootstrappedTrainingSamples[0][treeNode.numFeatures]
			return

		ifAlreadyClassified, treeNode.outputClass = self.checkIfAlreadyClassified(treeNode)
		if ifAlreadyClassified:
			print "Already classified: Case 2"
			treeNode.isLeaf = True
			return

		treeNode.featureSubsetIndices = self.selectFeatureSubsetRandomly()
		# DEBUG
#		treeNode.featureSubsetIndices = [100, 103, 161, 163, 197, 355, 460, 463, 507, 509]

		bestFeatureToSplitOn, bestFeatureSplitValue, bestFeatureSplitEntropy =  self.getBestFeatureToSplitOn(treeNode)
		# TODO: create a bound check, if you have the split on the extreme indices, that means already classified
		# so return immediately


		treeNode.splitFeatureIndex = bestFeatureToSplitOn
		treeNode.splitFeatureValue = bestFeatureSplitValue
		treeNode.splitFeatureEntropy = bestFeatureSplitEntropy

		print bestFeatureToSplitOn, bestFeatureSplitValue, bestFeatureSplitEntropy
		leftChildSamples, rightChildSamples = self.getSplitPopulation(treeNode)
		print "leftChildSamples:", leftChildSamples
		print "rightChildSamples:", rightChildSamples

		leftChildNode = TreeNode(leftChildSamples, self.numFeatures, len(leftChildSamples), self.numOutputClasses, treeNode.generation + 1)
		rightChildNode = TreeNode(rightChildSamples, self.numFeatures, len(rightChildSamples), self.numOutputClasses, treeNode.generation + 1)

		treeNode.leftChildNode = leftChildNode
		treeNode.rightChildNode = rightChildNode
		self.splitRecursively(leftChildNode)
		self.splitRecursively(rightChildNode)

		# given a split point, gives the user two different sets of data
	def getSplitPopulation(self, node):
		leftChildSamples, rightChildSamples = [], []
		globalIndex = node.splitFeatureIndex
		for x in node.bootstrappedTrainingSamples:
			if x[globalIndex] < node.splitFeatureValue: leftChildSamples.append(x)
			else: rightChildSamples.append(x)
		return leftChildSamples, rightChildSamples

	# given the feature indices, selects the best feature index to split on
	def getBestFeatureToSplitOn(self, node):
		bootstrappedFeatureVectors = node.bootstrappedFeatureVectors
		bootstrappedOutputVector = node.bootstrappedOutputVector
		featureSubsetIndices = node.featureSubsetIndices

		print "getBestFeatureToSplitOn()"
		for x in bootstrappedFeatureVectors: print x
		print "featureSubsetIndices: ", featureSubsetIndices
		featureSubsetEntropies = []
		featureSubsetSplitValues = []
		for i in range(0, len(featureSubsetIndices)):
			tryIndex = featureSubsetIndices[i]
			print "tryIndex", tryIndex
			featureMinEntropy, featureSplitValue = self.getMinEntropyOfFeature(bootstrappedFeatureVectors[tryIndex], bootstrappedOutputVector, self.numOutputClasses)
			featureSubsetEntropies.append(featureMinEntropy)
			featureSubsetSplitValues.append(featureSplitValue)

		print "featureSubsetEntropies:", featureSubsetEntropies
		print "featureSubsetSplitValues", featureSubsetSplitValues

		featureMinEntropy = min(featureSubsetEntropies)
		bestFeatureToSplitOn = featureSubsetEntropies.index(featureMinEntropy)
		bestFeatureSplitValue = featureSubsetSplitValues[bestFeatureToSplitOn]

#		print 'Best Split is possible with global feature:', featureSubsetIndices[bestFeatureToSplitOn]
#		print 'Which has an entropy of:', featureMinEntropy
#		print 'Best Split on this feature with value:', bestFeatureSplitValue
		return featureSubsetIndices[bestFeatureToSplitOn], bestFeatureSplitValue, featureMinEntropy

	def getMinEntropyOfFeature(self, featureVector, outputVector, numOutputClasses):
		print "getMinEntropyOfFeature()"
		# create feature vs output tuple
		featureOutputPair = [[featureVector[x], outputVector[x]] for x in range(0, len(featureVector))]
		featureOutputPair = sorted(featureOutputPair, key = lambda x: x[0])
		print "featureOutputPair", featureOutputPair

		#	for x in featureOutputPair: print x

		#	find the splitPoints that has the values where the changes has occurred,
		#	so these 'splitPoints' are best places to split
		searchOutput = featureOutputPair[0][1]
		splitPoints = []
		for i in range(0, len(featureVector)):
			if featureOutputPair[i][1] != searchOutput:
				splitPoints.append(i)
				searchOutput = featureOutputPair[i][1]

		#	now we have the splitPoints, so we need to find which of them gives
		#	us the highest entropy gain
		print "splitPoints:", splitPoints
		minEntropy, bestSplitIndex = self.getBestSplitAndMinEntropy(featureOutputPair, splitPoints, numOutputClasses)
		featureSplitValue = featureOutputPair[splitPoints[bestSplitIndex]][0]
#		print minEntropy, bestSplitIndex, splitPoints[bestSplitIndex], featureOutputPair[splitPoints[bestSplitIndex]][0]
		return minEntropy, featureSplitValue

	# calculate all the possible splits for a feature vector and then return the value of the best split
	def getBestSplitAndMinEntropy(self, featureOutputPair, splitPoints, numOutputClasses):
		print "getBestSplitAndMinEntropy()"
		numSamples = len(featureOutputPair)
#		print "numSamples:", numSamples

		entropies = []
		for index in splitPoints:
			# print "for the cut in index", index
			valueAtSplitPoint = featureOutputPair[index][0]
			numLessThanValueAtSplitPoint = 0
			numGreaterThanValueAtSplitPoint = 0
			for record in featureOutputPair:
				if record[0] <= valueAtSplitPoint: numLessThanValueAtSplitPoint += 1
				else: numGreaterThanValueAtSplitPoint += 1
			# print 'numLessThanValueAtSplitPoint:', numLessThanValueAtSplitPoint, 'numGreaterThanValueAtSplitPoint:', numGreaterThanValueAtSplitPoint
			# call for upper portion
			upperEntropyOfSplit = self.calcSplitEntropy(featureOutputPair, index, numOutputClasses, True)
			# call for lower portion
			lowerEntropyOfSplit = self.calcSplitEntropy(featureOutputPair, index, numOutputClasses, False)
			print "upperEntropyOfSplit:", upperEntropyOfSplit, "lowerEntropyOfSplit:", lowerEntropyOfSplit
			print "numLessThanValueAtSplitPoint:", numLessThanValueAtSplitPoint, "numGreaterThanValueAtSplitPoint:", numGreaterThanValueAtSplitPoint
			totalEntropy = (numLessThanValueAtSplitPoint * upperEntropyOfSplit + numGreaterThanValueAtSplitPoint * lowerEntropyOfSplit) / numSamples
#			print "totalEntropy", totalEntropy
			entropies.append(totalEntropy)

		print "entropies:", entropies
		minEntropy = min(entropies)
		minEntropyIndex = entropies.index(minEntropy)
		return minEntropy, minEntropyIndex

	# calculate entropy for each of the splits
	def calcSplitEntropy(self, featureOutputPairs, splitIndex, numOutputClasses, isUpperSplit):
#		print "calcSplitEntropy()"
#		print featureOutputPairs
		classCounts = [0 for i in range(0, numOutputClasses)]
		if isUpperSplit:
			for i in range(0, splitIndex): classCounts[featureOutputPairs[i][1]] += 1
		else:
			for i in range(splitIndex, len(featureOutputPairs)): classCounts[featureOutputPairs[i][1]] += 1
		totalClassCounts = sum(classCounts)

		splitEntropy = 0
		for x in classCounts:
			# if x == 0 we are in trouble, cause the log function will return infinity,
			# but we want actually zero so do a continue and skip
			if not x: continue
			probability = x / totalClassCounts
			splitEntropy += -(probability * log(probability, 2))

		return splitEntropy

	def checkIfAlreadyClassified(self, treeNode):
		print "checkIfAlreadyClassified()"
		outPutClasses = []
		for x in treeNode.bootstrappedTrainingSamples:
			if x[treeNode.numFeatures] not in outPutClasses: outPutClasses.append(x[treeNode.numFeatures])
		if len(outPutClasses) < 2: return True, outPutClasses[0]
		else: return False, None


class TreeNode:
	def __init__(self, bootstrappedTrainingSamples, numFeatures, numSamples, numOutputClasses, generation):
		self.numFeatures = numFeatures
		self.numSamples =  numSamples
		self.numOutputClasses = numOutputClasses
		self.isLeaf = False
		self.outputClass = None
		self.bootstrappedTrainingSamples = bootstrappedTrainingSamples
		self.bootstrappedFeatureVectors = []
		self.generation = generation

		self.splitFeatureIndex = 0
		self.splitFeatureValue = 0
		self.splitFeatureEntropy = 0

		self.calcFeatureVectors()
		self.bootstrappedOutputVector = [bootstrappedTrainingSamples[x][self.numFeatures] for x in range(0, self.numSamples)]
		self.featureSubsetIndices = []

		self.leftChildNode = None
		self.rightChildNode = None

	def calcFeatureVectors(self):
		print "calcFeatureVectors()"
		self.bootstrappedFeatureVectors = [[] for x in range (0, self.numFeatures)]
		for i in range(0, self.numSamples):
			for j in range(0, self.numFeatures):
				self.bootstrappedFeatureVectors[j].append(self.bootstrappedTrainingSamples[i][j])


''' The main algorithm for Regularized Random Forest, it creates a number of decision trees and then aggregates them
	to find the solution '''
class RegularizedRandomForest:
	def __init__(self, dataSet, numDecisionTrees):
		self.decisionTrees = []
		self.dataSet = dataSet
		self.numDecisionTrees = numDecisionTrees

	def populateDecisionTrees(self):
		for i in range(0, self.numDecisionTrees):
#			self.decisionTrees.append(createDecisionTree(self.dataSet))
			decisionTree = DecisionTree(dataSet)
			self.decisionTrees.append(decisionTree)

		# TODO do the usual stuffs (aggregation) with the decisionTrees

''' This class reads the file and crates a data matrix for further processing '''
class FileReader:

	def __init__(self, sharedFilePath, designFilePath):
		self.sharedFileData = []
		self.designFileData = []
		with open(sharedFilePath) as sharedFile:
			for line in sharedFile: self.sharedFileData.append(line.strip().split())
		with open(designFilePath) as designFile:
			for line in designFile: self.designFileData.append(line.strip().split())

	def getDataSet(self):
		self.dataSet = []
		self.dataSetClasses = {}
		for line in self.sharedFileData[1:]:
			self.dataSet.append([int(x) for x in line[3:]])
		i, j = 0, 0
		for line in self.designFileData:
			if line[1] not in self.dataSetClasses.keys():
				self.dataSetClasses[line[1]] = j
				j += 1
			self.dataSet[i].append(int(self.dataSetClasses[line[1]]))
			i += 1
		return self.dataSet

def readFileContents(fileName):
	data = []
	file = open(fileName)
	for x in file: data.append([int(y) for y in x.strip().split()])
	return data


if __name__ == "__main__":
	numDecisionTrees = 1

	# small-alter.txt has a modified dataset
#	dataSet = readFileContents('Datasets/small-alter.txt')

	mouseData = ["Datasets/final.an.0.03.subsample.0.03.pick.shared", "Datasets/mouse.sex_time.design"]
	sharedFilePath, designFilePath = mouseData
	fileReader = FileReader(sharedFilePath, designFilePath)
	dataSet = fileReader.getDataSet()
#	for x in dataSet: print dataSet

	regularizedRandomForest = RegularizedRandomForest(dataSet, numDecisionTrees)
	regularizedRandomForest.populateDecisionTrees()

#	createDecisionTree(dataSet)