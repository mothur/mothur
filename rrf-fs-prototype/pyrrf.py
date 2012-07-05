#!/usr/bin/python
from __future__ import division
from math import log, ceil, sqrt
import random

DEBUG_MODE = True

class AbstractDecisionTree(object):
	def __init__(self, baseDataSet, globalDiscardedFeatureIndices, optimumFeatureSubsetSelector):

		self.baseDataSet = baseDataSet
		self.numSamples = len(baseDataSet)
		self.numFeatures = len(baseDataSet[0]) - 1
		self.numOutputClasses = 0
		self.classes = []

		self.bootstrappedTrainingSamples = []
		self.bootstrappedTrainingSampleIndices = []
		self.bootstrappedTestSamples = []
		self.bootstrappedTestSampleIndices = []

		self.rootNode = None
		self.globalDiscardedFeatureIndices = globalDiscardedFeatureIndices

		self.optimumFeatureSubsetSize = optimumFeatureSubsetSelector.getOptimumFeatureSubsetSize(self.numFeatures)

		for i in range(0, self.numSamples):
			outcome = baseDataSet[i][-1]
			if outcome not in self.classes:
				self.classes.append(outcome)
				self.numOutputClasses += 1

	# this is the code which creates bootstrapped samples from given data
	def createBootStrappedSamples(self):
		isInTrainingSamples = [False for i in range(0, self.numSamples)]

		for i in range(0, self.numSamples):
			randomIndex = random.randint(0, self.numSamples - 1)
			self.bootstrappedTrainingSamples.append(self.baseDataSet[randomIndex])
			isInTrainingSamples[randomIndex] = True

		for i in range(0, self.numSamples):
			if isInTrainingSamples[i]:
				self.bootstrappedTrainingSampleIndices.append(i)
			else:
				self.bootstrappedTestSamples.append(self.baseDataSet[i])
				self.bootstrappedTestSampleIndices.append(i)

			#		for x in self.bootstrappedTrainingSamples: print x
			#		print "###"
		if DEBUG_MODE:
			for x in self.bootstrappedTestSamples: print x
		if DEBUG_MODE: print "self.bootstrappedTrainingSampleIndices:", self.bootstrappedTrainingSampleIndices
		if DEBUG_MODE: print "self.bootstrappedTestSampleIndices:", self.bootstrappedTestSampleIndices

	def getMinEntropyOfFeature(self, featureVector, outputVector, numOutputClasses):
		if DEBUG_MODE: print "getMinEntropyOfFeature()"
		# create feature vs output tuple
		featureOutputPair = zip(featureVector, outputVector)
		featureOutputPair = sorted(featureOutputPair, key = lambda x: x[0])
		if DEBUG_MODE: print "featureOutputPair", featureOutputPair

		#	for x in featureOutputPair: print x

		#	find the splitPoints that has the values where the changes has occurred,
		#	so these 'splitPoints' are best places to split
		#		searchOutput = featureOutputPair[0][1]
		splitPoints = []
		uniqueFeatureValues = [ featureOutputPair[0][0] ]

		# trying out a different implementation of split points calculation logic
		# this seems to be working fine
		# need to come up with a logic to speed this up
		for index, pair in enumerate(featureOutputPair):
			if pair[0] not in uniqueFeatureValues:
				uniqueFeatureValues.append(pair[0])
				splitPoints.append(index)

		#	now we have the splitPoints, so we need to find which of them gives
		#	us the highest entropy gain
		if DEBUG_MODE: print "splitPoints:", splitPoints

		# if we had too many zero containing split points, that means we ignored them all, and we have
		# now an empty split point list, so we need to check this before going further
		if not len(splitPoints):
			minEntropy = float("inf")
			bestSplitIndex = -1
			featureSplitValue = 0
		else:
			minEntropy, bestSplitIndex = self.getBestSplitAndMinEntropy(featureOutputPair, splitPoints, numOutputClasses)
			featureSplitValue = featureOutputPair[splitPoints[bestSplitIndex]][0]
		#		print minEntropy, bestSplitIndex, splitPoints[bestSplitIndex], featureOutputPair[splitPoints[bestSplitIndex]][0]

		return minEntropy, featureSplitValue


	# calculate all the possible splits for a feature vector and then return the value of the best split
	def getBestSplitAndMinEntropy(self, featureOutputPair, splitPoints, numOutputClasses):
		if DEBUG_MODE: print "getBestSplitAndMinEntropy()"
		numSamples = len(featureOutputPair)
		#		print "numSamples:", numSamples

		entropies = []
		for index in splitPoints:
			# print "for the cut in index", index
			valueAtSplitPoint = featureOutputPair[index][0]
			numLessThanValueAtSplitPoint = 0
			numGreaterThanValueAtSplitPoint = 0
			for record in featureOutputPair:
				if record[0] < valueAtSplitPoint: numLessThanValueAtSplitPoint += 1
				else: numGreaterThanValueAtSplitPoint += 1
			# print 'numLessThanValueAtSplitPoint:', numLessThanValueAtSplitPoint, 'numGreaterThanValueAtSplitPoint:', numGreaterThanValueAtSplitPoint
			# call for upper portion
			upperEntropyOfSplit = self.calcSplitEntropy(featureOutputPair, index, numOutputClasses, True)
			# call for lower portion
			lowerEntropyOfSplit = self.calcSplitEntropy(featureOutputPair, index, numOutputClasses, False)
			if DEBUG_MODE: print "upperEntropyOfSplit:", upperEntropyOfSplit, "lowerEntropyOfSplit:", lowerEntropyOfSplit
			if DEBUG_MODE: print "numLessThanValueAtSplitPoint:", numLessThanValueAtSplitPoint, "numGreaterThanValueAtSplitPoint:", numGreaterThanValueAtSplitPoint
			totalEntropy = (numLessThanValueAtSplitPoint * upperEntropyOfSplit + numGreaterThanValueAtSplitPoint * lowerEntropyOfSplit) / numSamples
			#			print "totalEntropy", totalEntropy
			entropies.append(totalEntropy)

		if DEBUG_MODE: print "entropies:", entropies
		minEntropy = min(entropies)
		minEntropyIndex = entropies.index(minEntropy)
		return minEntropy, minEntropyIndex


# This is the regularized version of the DecisionTree class
class RegularizedDecisionTree(AbstractDecisionTree):
	def __init__(self, baseDataSet, globalDiscardedFeatureIndices, optimumFeatureSubsetSelector):

		# calling to super-class's ctor
		super(RegularizedDecisionTree, self).__init__(baseDataSet, globalDiscardedFeatureIndices, optimumFeatureSubsetSelector)

		self.createBootStrappedSamples()
		self.buildDecisionTree()

	def buildDecisionTree(self):
		if DEBUG_MODE: print "buildDecisionTree()"
		treeNode = TreeNode(self.bootstrappedTrainingSamples, self.globalDiscardedFeatureIndices, self.numFeatures, self.numSamples, self.numOutputClasses, 0)
		self.rootNode = treeNode

		penalty = 1
		featureSubset = []
		self.splitRecursively(treeNode, featureSubset, penalty)

		# TODO: implement the printTree method
#		if DEBUG_MODE: self.printTree(treeNode, "root")

	def splitRecursively(self, treeNode, featureSubset, penalty):
		if DEBUG_MODE: print "splitRecursively()"

		bestGain = 0
		numberOfNewFeaturesTested = 0
		penalizedGains = [0 for x in range(0, self.numFeatures)]

		# TODO: how about we shuffle the value of the indices rather than using increasing
		# value? the randomization might give it a quality boost
		for x in range(0, self.numFeatures):
			if x in featureSubset:
				# calculate the gain by using x
#				node.bootstrappedFeatureVectors
				pass
			else:
				pass


#		# return immediately if this node is just a leaf, no recursion is needed
#		if treeNode.numSamples < 2:
#			if DEBUG_MODE: print "Already classified: Case 1"
#			treeNode.isLeaf = True
#			treeNode.outputClass = treeNode.bootstrappedTrainingSamples[0][treeNode.numFeatures]
#			return
#
#		ifAlreadyClassified, treeNode.outputClass = self.checkIfAlreadyClassified(treeNode)
#		if ifAlreadyClassified:
#			if DEBUG_MODE: print "Already classified: Case 2"
#			treeNode.isLeaf = True
#			return
#
#		#		nullFeatureIndices = self.findNullFeatures(treeNode.bootstrappedTrainingSamples)
#		treeNode.featureSubsetIndices = self.selectFeatureSubsetRandomly(self.globalDiscardedFeatureIndices, treeNode.localDiscardedFeatureIndices)
#		if DEBUG_MODE: print "discardedFeatureIndices:", globalDiscardedFeatureIndices
#		if DEBUG_MODE: print "featureSubsetIndices:", treeNode.featureSubsetIndices
#		# DEBUG
#		#		treeNode.featureSubsetIndices = [100, 103, 161, 163, 197, 355, 460, 463, 507, 509]
#
#		bestFeatureToSplitOn, bestFeatureSplitValue, bestFeatureSplitEntropy =  self.getBestFeatureToSplitOn(treeNode)
#		# TODO: create a bound check, if you have the split on the extreme indices, that means already classified
#		# so return immediately
#
#
#		treeNode.splitFeatureIndex = bestFeatureToSplitOn
#		treeNode.splitFeatureValue = bestFeatureSplitValue
#		treeNode.splitFeatureEntropy = bestFeatureSplitEntropy
#
#		if DEBUG_MODE: print "bestFeatureToSplitOn:", bestFeatureToSplitOn, "bestFeatureSplitValue:", bestFeatureSplitValue, "bestFeatureSplitEntropy:", bestFeatureSplitEntropy
#		leftChildSamples, rightChildSamples = self.getSplitPopulation(treeNode)
#		# print "leftChildSamples:", leftChildSamples
#		# print "rightChildSamples:", rightChildSamples
#
#		leftChildNode = TreeNode(leftChildSamples, self.globalDiscardedFeatureIndices, self.numFeatures, len(leftChildSamples), self.numOutputClasses, treeNode.generation + 1)
#		rightChildNode = TreeNode(rightChildSamples, self.globalDiscardedFeatureIndices, self.numFeatures, len(rightChildSamples), self.numOutputClasses, treeNode.generation + 1)
#
#		treeNode.leftChildNode = leftChildNode
#		treeNode.rightChildNode = rightChildNode
#		self.splitRecursively(leftChildNode)
#		self.splitRecursively(rightChildNode)




# The main algorithm resides here, it the manipulator class
class DecisionTree(AbstractDecisionTree):
	def __init__(self, baseDataSet, globalDiscardedFeatureIndices, optimumFeatureSubsetSelector):

		# calling to super-class's ctor
		super(DecisionTree, self).__init__(baseDataSet, globalDiscardedFeatureIndices, optimumFeatureSubsetSelector)

		self.variableImportanceList = [0 for x in range(0, self.numFeatures)]
		self.outOfBagEstimates = {}

		self.createBootStrappedSamples()
		self.buildDecisionTree()


	def findNullFeatures(self, bootstrappedTrainingSamples):
		print "findNullFeatures()"
		# we can easily get the featureVectors by transposing, this zip functions is an easy
		# shortcut for transposing a matrix
		featureVectors = zip(*bootstrappedTrainingSamples)
		nullFeatureIndices = []
		for index, value in enumerate(featureVectors):
			total = sum(value)
			if not total: nullFeatureIndices.append(index)
		return nullFeatureIndices


	# randomly select log2(totalFeatures) number features for each node
	def selectFeatureSubsetRandomly(self, globalDiscardedFeatureIndices, localDiscardedFeatureIndices):
		if DEBUG_MODE: print "selectFeatureSubsetRandomly()"

		featureSubsetIndices = []

		combinedDiscardedFeatureIndices = []
		combinedDiscardedFeatureIndices.extend(globalDiscardedFeatureIndices)
		combinedDiscardedFeatureIndices.extend(localDiscardedFeatureIndices)

		numberOfRemainingSuitableFeatures = self.numFeatures - len(combinedDiscardedFeatureIndices)

		if numberOfRemainingSuitableFeatures < self.optimumFeatureSubsetSize: currentFeatureSubsetSize = numberOfRemainingSuitableFeatures
		else: currentFeatureSubsetSize = self.optimumFeatureSubsetSize

		while len(featureSubsetIndices) < currentFeatureSubsetSize:
			randomIndex = random.randint(0, self.numFeatures - 1)
			# TODO the loop goes infinite here since it cannot find remaining features, need a way to break the loop
			if (randomIndex not in featureSubsetIndices) and (randomIndex not in combinedDiscardedFeatureIndices):
				featureSubsetIndices.append(randomIndex)
#		print 'returning from selectFeatureSubsetRandomly()'
		return sorted(featureSubsetIndices)

	def buildDecisionTree(self):
		if DEBUG_MODE: print "buildDecisionTree()"
		treeNode = TreeNode(self.bootstrappedTrainingSamples, self.globalDiscardedFeatureIndices, self.numFeatures, self.numSamples, self.numOutputClasses, 0)
		self.rootNode = treeNode
		self.splitRecursively(treeNode)
		if DEBUG_MODE: self.printTree(treeNode, "root")

	def printTree(self, treeNode, caption):
		tabs = ""
		for i in range(0, treeNode.generation): tabs += '\t'

		if not treeNode.isLeaf:
			print tabs + caption + ' [ gen: ' + str(treeNode.generation) + ' ] ( ' + str(treeNode.splitFeatureValue) + ' < X'+ str(treeNode.splitFeatureIndex) + ' )'
			self.printTree(treeNode.leftChildNode, 'leftChild')
			self.printTree(treeNode.rightChildNode, 'rightChild')
		else:
			print tabs + caption + ' [ gen: ' + str(treeNode.generation) + ' ] ( classified to: ' + str(treeNode.outputClass) +', samples: ' + str(treeNode.numSamples) +' )'


	def splitRecursively(self, rootNode):
		if DEBUG_MODE: print "splitRecursively()"

		# return immediately if this node is just a leaf, no recursion is needed
		if rootNode.numSamples < 2:
			if DEBUG_MODE: print "Already classified: Case 1"
			rootNode.isLeaf = True
			rootNode.outputClass = rootNode.bootstrappedTrainingSamples[0][rootNode.numFeatures]
			return

		ifAlreadyClassified, rootNode.outputClass = self.checkIfAlreadyClassified(rootNode)
		if ifAlreadyClassified:
			if DEBUG_MODE: print "Already classified: Case 2"
			rootNode.isLeaf = True
			return
		
#		nullFeatureIndices = self.findNullFeatures(treeNode.bootstrappedTrainingSamples)
		rootNode.featureSubsetIndices = self.selectFeatureSubsetRandomly(self.globalDiscardedFeatureIndices, rootNode.localDiscardedFeatureIndices)
		if DEBUG_MODE: print "discardedFeatureIndices:", self.globalDiscardedFeatureIndices
		if DEBUG_MODE: print "featureSubsetIndices:", rootNode.featureSubsetIndices
		# DEBUG
#		treeNode.featureSubsetIndices = [100, 103, 161, 163, 197, 355, 460, 463, 507, 509]

		bestFeatureToSplitOn, bestFeatureSplitValue, bestFeatureSplitEntropy =  self.getBestFeatureToSplitOn(rootNode)
		# TODO: create a bound check, if you have the split on the extreme indices, that means already classified
		# so return immediately


		rootNode.splitFeatureIndex = bestFeatureToSplitOn
		rootNode.splitFeatureValue = bestFeatureSplitValue
		rootNode.splitFeatureEntropy = bestFeatureSplitEntropy

		if DEBUG_MODE: print "bestFeatureToSplitOn:", bestFeatureToSplitOn, "bestFeatureSplitValue:", bestFeatureSplitValue, "bestFeatureSplitEntropy:", bestFeatureSplitEntropy
		leftChildSamples, rightChildSamples = self.getSplitPopulation(rootNode)
		# print "leftChildSamples:", leftChildSamples
		# print "rightChildSamples:", rightChildSamples

		leftChildNode = TreeNode(leftChildSamples, self.globalDiscardedFeatureIndices, self.numFeatures, len(leftChildSamples), self.numOutputClasses, rootNode.generation + 1)
		rightChildNode = TreeNode(rightChildSamples, self.globalDiscardedFeatureIndices, self.numFeatures, len(rightChildSamples), self.numOutputClasses, rootNode.generation + 1)

		rootNode.leftChildNode = leftChildNode
		rootNode.leftChildNode.parentNode = rootNode
		rootNode.rightChildNode = rightChildNode
		rootNode.rightChildNode.parentNode = rootNode

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

		if DEBUG_MODE: print "getBestFeatureToSplitOn()"
#		for x in bootstrappedFeatureVectors: print x
#		print "featureSubsetIndices: ", featureSubsetIndices
		featureSubsetEntropies = []
		featureSubsetSplitValues = []
		for i in range(0, len(featureSubsetIndices)):
			tryIndex = featureSubsetIndices[i]
			if DEBUG_MODE: print "trying feature of index:", tryIndex
			featureMinEntropy, featureSplitValue = self.getMinEntropyOfFeature(bootstrappedFeatureVectors[tryIndex], bootstrappedOutputVector, self.numOutputClasses)
			featureSubsetEntropies.append(featureMinEntropy)
			featureSubsetSplitValues.append(featureSplitValue)

		if DEBUG_MODE: print "featureSubsetEntropies:", featureSubsetEntropies
		if DEBUG_MODE: print "featureSubsetSplitValues", featureSubsetSplitValues

		featureMinEntropy = min(featureSubsetEntropies)
		bestFeatureToSplitOn = featureSubsetEntropies.index(featureMinEntropy)
		bestFeatureSplitValue = featureSubsetSplitValues[bestFeatureToSplitOn]

#		print 'Best Split is possible with global feature:', featureSubsetIndices[bestFeatureToSplitOn]
#		print 'Which has an entropy of:', featureMinEntropy
#		print 'Best Split on this feature with value:', bestFeatureSplitValue
		return featureSubsetIndices[bestFeatureToSplitOn], bestFeatureSplitValue, featureMinEntropy



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
		if DEBUG_MODE: print "checkIfAlreadyClassified()"
		if DEBUG_MODE: print "len(bootstrappedTrainingSamples):", len(treeNode.bootstrappedTrainingSamples)
#		if len(treeNode.bootstrappedTrainingSamples) < 10:
#			print treeNode.bootstrappedTrainingSamples
		outPutClasses = []
		for i, x in enumerate(treeNode.bootstrappedTrainingSamples):
			# print "index:", i
			if x[treeNode.numFeatures] not in outPutClasses:
				if DEBUG_MODE: print "appending: ", x[treeNode.numFeatures]
				outPutClasses.append(x[treeNode.numFeatures])
		if len(outPutClasses) < 2: return True, outPutClasses[0]
		else: return False, None

	def calcTreeVariableImportanceAndError(self):
		if DEBUG_MODE: print "calcTreeVariableImportanceAndError()"
		numCorrect, treeErrorRate = self.calcTreeErrorRate()
		print "len(self.bootstrappedTestSamples):", len(self.bootstrappedTestSamples)
		print "numCorrect:", numCorrect
		print "treeErrorRate:", treeErrorRate


		for i in range(0, self.numFeatures):
			randomlySampledTestData = self.randomlyShuffleAttribute(self.bootstrappedTestSamples, i)

			numCorrectAfterShuffle = 0
			for shuffledSample in randomlySampledTestData:
				actualSampleOutputClass = shuffledSample[self.numFeatures]
				predictedSampleOutputClass = self.evaluateSample(shuffledSample)
				if actualSampleOutputClass == predictedSampleOutputClass:
					numCorrectAfterShuffle += 1
			self.variableImportanceList[i] += (numCorrect - numCorrectAfterShuffle)

#		print "self.variableImportanceList:", self.variableImportanceList

		variableRanks = []
		for i, x in enumerate(self.variableImportanceList):
			if x > 0: variableRanks.append([i, x])
		variableRanks = sorted(variableRanks, key = lambda x: x[1], reverse = True)
		print "variableRanks:", variableRanks


	# given a sample and a featureIndex, it randomly shuffles that only specified feature
	# and returns the whole dataset
	def randomlyShuffleAttribute(self, samples, featureIndex):
		featureVectors = [list(x) for x in zip(*samples)]
		random.shuffle(featureVectors[featureIndex])
		shuffledSample = zip(*featureVectors)
		return shuffledSample

	# uses the evaluateSample() function to calculate the error rate of the tree
	def calcTreeErrorRate(self):
		numCorrect = 0
		for i, testSample in enumerate(self.bootstrappedTestSamples):
			testSampleIndex = self.bootstrappedTestSampleIndices[i]

			actualSampleOutputClass = testSample[self.numFeatures]
			predictedSampleOutputClass = self.evaluateSample(testSample)
			if actualSampleOutputClass == predictedSampleOutputClass: numCorrect += 1

			self.outOfBagEstimates[testSampleIndex] = predictedSampleOutputClass

		treeErrorRate = 1 - (numCorrect / len(self.bootstrappedTestSamples))
		return numCorrect, treeErrorRate

	# this function will be used to get the predicted output by giving it a test data
	def evaluateSample(self, testSample):
		node = self.rootNode
		while True:
			if node.isLeaf: return node.outputClass

			sampleSplitFeatureValue = testSample[node.splitFeatureIndex]
			if sampleSplitFeatureValue < node.splitFeatureValue: node = node.leftChildNode
			else: node = node.rightChildNode

	def purgeDataSetsFromTree(self):
		self.purgeTreeNodesDataRecursively(self.rootNode)

	# once the calculation has been done, we do not need to have the samples associated with the TreeNodes, this
	# hogs a lot of data. So recursively delete the values associated with the TreeNode, just keep the values that
	# are important for classification
	def purgeTreeNodesDataRecursively(self, treeNode):
		treeNode.bootstrappedTrainingSamples = None
		treeNode.bootstrappedFeatureVectors = None
		treeNode.bootstrappedOutputVector = None
		treeNode.localDiscardedFeatureIndices = None
		treeNode.globalDiscardedFeatureIndices = None

		if treeNode.leftChildNode is not None: self.purgeTreeNodesDataRecursively(treeNode.leftChildNode)
		if treeNode.rightChildNode is not None: self.purgeTreeNodesDataRecursively(treeNode.rightChildNode)

class TreeNode(object):
	def __init__(self, bootstrappedTrainingSamples, globalDiscardedFeatureIndices, numFeatures, numSamples, numOutputClasses, generation):
		self.numFeatures = numFeatures
		self.numSamples =  numSamples
		self.numOutputClasses = numOutputClasses
		self.isLeaf = False
		self.outputClass = None
		self.bootstrappedTrainingSamples = bootstrappedTrainingSamples
		self.globalDiscardedFeatureIndices = globalDiscardedFeatureIndices
		self.generation = generation

		self.splitFeatureIndex = 0
		self.splitFeatureValue = 0
		self.splitFeatureEntropy = 0

		self.bootstrappedFeatureVectors = [list(x) for x in zip(*self.bootstrappedTrainingSamples)]
		self.bootstrappedOutputVector = [bootstrappedTrainingSamples[x][self.numFeatures] for x in range(0, self.numSamples)]
		self.featureSubsetIndices = []

		self.leftChildNode = None
		self.rightChildNode = None
		self.parentNode = None

		self.localDiscardedFeatureIndices = []

		# call some helper functions
		self.createLocalDiscardedFeatureList()

	def createLocalDiscardedFeatureList(self):
		if DEBUG_MODE: print "createLocalDiscardedFeatureList()"
		for i, x in enumerate(self.bootstrappedFeatureVectors):
			if i not in self.globalDiscardedFeatureIndices and Utils.getStandardDeviation(x) <= 0:
				self.localDiscardedFeatureIndices.append(i)
		if DEBUG_MODE: print self.localDiscardedFeatureIndices

class AbstractRandomForest(object):
	def __init__(self, dataSet, numDecisionTrees):
		self.decisionTrees = []
		self.dataSet = dataSet
		self.numDecisionTrees = numDecisionTrees

		self.globalDiscardedFeatureIndices = self.getGlobalDiscardedFeatureIndices()

		self.numSamples = len(dataSet)
		self.numFeatures = len(dataSet[0]) - 1

		self.globalOutOfBagEstimates = {}
		self.globalVariableImportanceList = [0 for x in range(0, self.numFeatures)]

	def getGlobalDiscardedFeatureIndices(self):
		featureVectors = zip(*self.dataSet)[:-1]
		globalDiscardedFeatureIndices = []
		for index, featureVector in enumerate(featureVectors):
			total = sum(featureVector)
			zeroCount = featureVector.count(0)
			if Utils.getStandardDeviation(featureVector) <= 0:
				globalDiscardedFeatureIndices.append(index)
		if DEBUG_MODE: print 'number of global discarded features:', len(globalDiscardedFeatureIndices)
		if DEBUG_MODE: print 'total features:', len(featureVectors)
		return globalDiscardedFeatureIndices

class RegularizedRandomForest(AbstractRandomForest):
	# TODO: finish the implementation
	def populateDecisionTrees(self):
		print "populateDecisionTrees()"
		for i in range(0, self.numDecisionTrees):
			print "Creating", i, "(th) Decision tree"
			decisionTree = RegularizedDecisionTree(dataSet, self.globalDiscardedFeatureIndices, OptimumFeatureSubsetSelector('squareRoot'))
#			decisionTree.calcTreeVariableImportanceAndError()
#
#			self.updateGlobalOutOfBagEstimates(decisionTree)
#
#			decisionTree.purgeDataSetsFromTree()
#
#			self.decisionTrees.append(decisionTree)
#
#		# TODO do the usual stuffs (aggregation) with the decisionTrees
#		if DEBUG_MODE: print "self.globalOutOfBagEstimates:", self.globalOutOfBagEstimates

		pass


''' The main algorithm for Regularized Random Forest, it creates a number of decision trees and then aggregates them
	to find the solution '''
class RandomForest(AbstractRandomForest):
	def updateGlobalOutOfBagEstimates(self, decisionTree):
		for indexOfSample, predictedOutcomeOfSample in decisionTree.outOfBagEstimates.iteritems():
			if not self.globalOutOfBagEstimates.has_key(indexOfSample):
				self.globalOutOfBagEstimates[indexOfSample] = [0 for i in range(0, decisionTree.numOutputClasses)]

			self.globalOutOfBagEstimates[indexOfSample][predictedOutcomeOfSample] += 1

	def calcForrestVariableImportance(self):
		print "calcForrestVariableImportance()"

		# make the sum of all the importance
		for decisionTree in self.decisionTrees:
			for i in range(0, self.numFeatures):
				self.globalVariableImportanceList[i] += decisionTree.variableImportanceList[i]

		# take the average of the trees
		for i in range(0, self.numFeatures):
			self.globalVariableImportanceList[i] /= self.numDecisionTrees

		globalVariableRanks = []
		for i, x in enumerate(self.globalVariableImportanceList):
			if x > 0: globalVariableRanks.append([i, x])
		globalVariableRanks = sorted(globalVariableRanks, key = lambda x: x[1], reverse = True)
		print "globalVariableRanks:", globalVariableRanks

	def calcForrestErrorRate(self):
		print "calcForrestErrorRate()"
		numCorrect = 0
		for indexOfSample, predictedOutComes in self.globalOutOfBagEstimates.iteritems():
			majorityVotedOutcome = predictedOutComes.index(max(predictedOutComes))
			realOutCome = self.dataSet[indexOfSample][self.numFeatures]

			if majorityVotedOutcome == realOutCome: numCorrect += 1
#			print "realOutCome", realOutCome, "majorityVotedOutcome", majorityVotedOutcome
		print "len(self.globalOutOfBagEstimates):", len(self.globalOutOfBagEstimates)
		print "numCorrect", numCorrect

		forrestErrorRate = 1 - (numCorrect / len(self.globalOutOfBagEstimates))
		print "forrestErrorRate:", forrestErrorRate

	def populateDecisionTrees(self):
		print "populateDecisionTrees()"
		for i in range(0, self.numDecisionTrees):
			print "Creating", i, "(th) Decision tree"
			decisionTree = DecisionTree(dataSet, self.globalDiscardedFeatureIndices, OptimumFeatureSubsetSelector('log2'))
			decisionTree.calcTreeVariableImportanceAndError()

			self.updateGlobalOutOfBagEstimates(decisionTree)

			decisionTree.purgeDataSetsFromTree()

			self.decisionTrees.append(decisionTree)

		# TODO do the usual stuffs (aggregation) with the decisionTrees
		if DEBUG_MODE: print "self.globalOutOfBagEstimates:", self.globalOutOfBagEstimates


class OptimumFeatureSubsetSelector(object):
	def __init__(self, selectionType = 'log2'):
		self.selectionType = selectionType

	def getOptimumFeatureSubsetSize(self, numFeatures):
		if DEBUG_MODE: print 'getOptimumFeatureSubsetSize()'
		if self.selectionType == 'log2':
			if DEBUG_MODE: print 'Using log2 mode'
			return int(ceil(log(numFeatures, 2)))
		elif self.selectionType == 'squareRoot':
			if DEBUG_MODE: print 'Using squareRoot mode'
			return int(ceil(sqrt(numFeatures)))
		else: return None


class FileReaderFactory(object):
	def __init__(self, fileType = "matrix", sharedFilePath = None, designFilePath = None, matrixFilePath = None):
		self.sharedFilePath = sharedFilePath
		self.designFilePath = designFilePath
		self.matrixFilePath = matrixFilePath

		self.fileReader = None

		if fileType == "matrix":
			self.fileReader = MatrixFileReader(self.matrixFilePath)
		elif fileType == "sharedAndDesign":
			self.fileReader = SharedAndDesignFileReader(self.sharedFilePath, designFilePath)

	def getFileReader(self):
		return self.fileReader

class MatrixFileReader(object):
	def __init__(self, matrixFilePath):
		self.matrixFilePath = matrixFilePath

	def getDataSetFromFile(self):
		data = []
		file = open(self.matrixFilePath)
		for line in file: data.append([int(y) for y in line.strip().split()])
		return data


''' This class reads the file and crates a data matrix for further processing '''
class SharedAndDesignFileReader(object):

	def __init__(self, sharedFilePath, designFilePath):
		self.sharedFileData = []
		self.designFileData = []
		with open(sharedFilePath) as sharedFile:
			for line in sharedFile: self.sharedFileData.append(line.strip().split())
		with open(designFilePath) as designFile:
			for line in designFile: self.designFileData.append(line.strip().split())

	def getDataSetFromFile(self):
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

class Utils(object):
	# standard deviation calculation function
	@staticmethod
	def getStandardDeviation(featureVector):
		n = len(featureVector)
		if not n:
			# standard deviation cannot be negative, this special value is returned to let the caller
			# function know that the list is empty
			return -1
		avg = sum(featureVector) / n
		standardDeviation = sqrt(sum([ (x - avg) ** 2 for x in featureVector ]) / n)
		return standardDeviation


if __name__ == "__main__":
	numDecisionTrees = 1

	# example of matrix file reading
#	fileReaderFactory = fileReaderFactory(fileType = 'matrix', matrixFilePath = 'Datasets/small-alter.txt');
	fileReaderFactory = FileReaderFactory(fileType = 'matrix', matrixFilePath = 'Datasets/inpatient.final.an.0.03.subsample.avg.matrix');
#	fileReaderFactory = FileReaderFactory(fileType = 'matrix', matrixFilePath = 'Datasets/outin.final.an.0.03.subsample.avg.matrix');

	# example of shared and design file reading
#	fileReaderFactory = FileReaderFactory(fileType='sharedAndDesign', sharedFilePath='Datasets/final.an.0.03.subsample.0.03.pick.shared', designFilePath='Datasets/mouse.sex_time.design')

	dataSet = fileReaderFactory.getFileReader().getDataSetFromFile()

	randomForest = RandomForest(dataSet, numDecisionTrees)
	randomForest.populateDecisionTrees()
	randomForest.calcForrestErrorRate()
	randomForest.calcForrestVariableImportance()

#	regularizedRandomForest = RegularizedRandomForest(dataSet, numDecisionTrees)
#	regularizedRandomForest.populateDecisionTrees()

#	createDecisionTree(dataSet)