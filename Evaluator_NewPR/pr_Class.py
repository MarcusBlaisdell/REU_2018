
import random
import time
from functions import sign

class pr_Class():
    ### global variables:

    T = 10 # default value of T
    k = 1 # default value of k
    kmer = 11 # default k-mer size is 11
    w = {} # weight is a dictionary
    eta = 1 # eta is learning rate
    trainMistakes = 0
    trainTotal = 0
    testMistakes = 0
    testTotal = 0
    trainPrecision = 0.0
    trainRecall = 0.0
    trainnpr = 0 # numerator for trainPrecision and trainRecall
    traindp = 0 # denominator for trainPrecision
    trainF1 = 0.0
    testPrecision = 0.0
    testRecall = 0.0
    testnpr = 0 # numerator for precision and recall
    testdp = 0 # denominator for Precision
    testF1 = 0.0
    trainGood = 0
    testGood = 0
    trainAccuracy = 0.0
    testAccuracy = 0.0
    tp = 0
    fp = 0
    tn = 0
    fn = 0
    scores = []
    labels = []
    target_recall = 0.98

    def __init__(self):
        pass

    ### perceptron function:

    def perceptron(self, trainDataSet):
        self.trainTotal += len(trainDataSet)

        for l in range (len(trainDataSet)):
            xit = trainDataSet[l][0]
            yit = trainDataSet[l][1]
            if yit == 1:
                self.trainGood += 1 ### increment count of actual good

            ### make the prediction: yHat = y*(<w,x>)
            ### changed to yHat = <w,x> to make the precision/recall calculations easier

            ### excluding the yit product lets me evaluate the relevant/irrelevant
            ### values without decoding them from a generic good prediction vs. bad prediction
            ### if yHat is a positive value, then it has predicted that the kmer is relevant
            ### if it is negative, it has predicted that the kmer is irrelevant

            theScore = self.dotProd(xit)
            self.scores.append(theScore)
            if yit == 1:
                self.labels.append(yit)
            else:
                self.labels.append(0)

            yHat = sign(theScore)

            #yHat = self.dotProd(xit) # Method U, real-valued yHat

            if yHat > 0: # if we predicted good:
                self.traindp += 1 # +1 prediction count
                if yit == 1: # if it is actually good:
                    self.trainnpr += 1 # correct +1 prediction
                else: # if it is actually false:
                    self.trainMistakes += 1
                    self.updateWeight(xit, yit, yHat)
            else: # if we predicted bad:
                if yit == 1: # if it is actually good:
                    self.trainMistakes += 1
                    self.updateWeight(xit, yit, yHat)

    ### end perceptron function

    ### use new weight to test accuracy on test dataList

    def testWeight (self, testDataSet, t):

        self.testTotal += len(testDataSet)

        ### evaluate all test samples:

        for i in range(len(testDataSet)):
            xit = testDataSet[i][0]
            yit = testDataSet[i][1]
            if yit == 1:
                self.testGood += 1 ### increment count of actual good

            yHat = sign(self.dotProd(xit)) # Method D

            ### if the value is actually good,
            ### and we predicted good, increment npr
            ### which is the # lines predicted as good that are actually good
            if yHat > 0: # if we predict good:
                self.testdp += 1 # +1 prediction count
                if yit == 1: #if it is actually good:
                    self.testnpr += 1 # correct +1 prediction
                    self.tp += 1 # increment true positive
                else: # if it is actually bad:
                    self.testMistakes += 1
                    self.fp += 1 # increment false positive

            else: ## if we predict bad:
                if yit == 1: #if it is actually good:
                    self.testMistakes += 1
                    self.fn += 1 # increment false FalseNegative
                else: #if it is actually bad:
                    self.tn += 1 # increment true negative
    ### end testWeight function

    ### function dotProd(), dot product by index
    def dotProd(self, xArray):
        #startTime = time.time()
        result = 0.0

        for element in xArray:
            if self.w.get(element[0], '--') == '--':
                self.w[element[0]] = 0
            else:
                result += self.w[element[0]] * element[1]

        #endTime = time.time()
        #print 'dotProd2: ', endTime - startTime

        return result

    ### end dotProd() function

    ### updateWeight function:
    def updateWeight(self, xarray, yStar, yHat):
        #startTime = time.time()

        ### update the weights in the weight vector that are in xit
        for element in xarray:
            ### if it is positive, and we predicted negative, we need to increase the weights
            ### minus a minus will do that (+1) - (-x) = 1 + x
            ### if it is negative and we predicted positive, we need to decrease the weights
            ### minus a positive will do that (-1) - (+x) = -1 -x

            #newY = yStar - yHat
            #self.w[element[0]] += self.eta * newY * element[1]
            self.w[element[0]] += self.eta * yStar * element[1]

    ### end updateWeight function
