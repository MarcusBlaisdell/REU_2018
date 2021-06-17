##########################################
###  Test
### class
### Marcus Blaisdell
##########################################

import random
import time
from functions import sign

class pr_Class():
    ### global variables:
    trainData = []
    validationData = []
    testData = []
    T = 1 # default value of T
    k = 1 # default value of k
    kmer = 11 # default k-mer size is 11
    #Tlist = [5, 10, 100, 1000, 5000, 10000]
    Tlist = [20]
    #klist = [1, 5, 10]
    klist = [1]
    w = {} # weight is a dictionary
    eta = 1.0 # eta is learning rate
    trainMistakes = 0
    trainTotal = 0
    testMistakes = 0
    testTotal = 0
    validationMistakes = 0
    validationTotal = 0
    trainPrecision = 0.0
    trainRecall = 0.0
    trainnpr = 0 # numerator for trainPrecision and trainRecall
    traindp = 0 # denominator for trainPrecision
    trainF1 = 0.0
    validationPrecision = 0.0
    validationRecall = 0.0
    validationnpr = 0 # numerator for precision and recall
    validationdp = 0 # denominator for Precision
    validationF1 = 0.0
    testPrecision = 0.0
    testRecall = 0.0
    testnpr = 0 # numerator for precision and recall
    testdp = 0 # denominator for Precision
    testF1 = 0.0
    trainGood = 0
    validationGood = 0
    testGood = 0
    trainAccuracy = 0.0
    validationAccuracy = 0.0
    testAccuracy = 0.0

    def __init__(self):
        pass

    ### perceptron function:

    #def perceptron(self, trainDataSet, testDataSet, outFile):
    def perceptron(self, trainDataSet):
        #self.trainTotal = len(trainDataSet)

        xit = [] # x-sub-i-sub-t, the training vector for the current iterations
        yit = 0 # y-sub-i-sub-t, the label for the current training vector (yStar)
        yHat = 0

        for l in range (len(trainDataSet)):
            self.trainTotal += 1
            xit = trainDataSet[l][0]
            yit = trainDataSet[l][1]

            ### make the prediction: yHat = y*(<w,x>)
            #yHat = yit * self.dotProd(xit) # Method U
            #yHat = yit * (self.dotProd(xit) + self.b) # Method V
            #yHat = sign(yit * (self.dotProd(xit) + self.b)) # Method W
            yHat = sign(yit * self.dotProd(xit)) # Method X

            #print 'initial eval: ', yit, ' : ', yHat

            ### update weight accordingly,
            ### if predicted value and actual value don't match,
            ## update weight,
            ### if they do match no update required
            if yit == 1:
                self.trainGood += 1
                if yHat > 0:
                    self.trainnpr += 1
                    self.traindp += 1
                else:
                    #print yit, ' : ', yHat
                    self.updateWeight(xit, yit)
                    self.trainMistakes += 1
            else:
                if yHat > 0:
                    #print yit, ' : ', yHat
                    self.trainMistakes += 1
                    self.updateWeight(xit, yit)
                else:
                    self.traindp += 1

            ### end train loop, trains on each sample in trainData

        ### end iteration loop

    ### end perceptron function

    ### use new weight to test accuracy on test dataList

    #def testWeight (self, trainDataSet, testDataSet, outFile, t):
    def testWeight (self, testDataSet, t):
        # reset mistakes count so each iteration starts at 0

        xit = [] # x-sub-i-sub-t, the training vector for the current iterations
        yit = 0 # y-sub-i-sub-t, the label for the current training vector (yStar)
        yHat = 0

        ### evaluate all test samples:

        for i in range(len(testDataSet)):
            self.testTotal += 1
            xit = testDataSet[i][0]
            yit = testDataSet[i][1]

            ### make the prediction:
            #yHat = yit * self.dotProd(xit) # Method A
            #yHat = yit * (self.dotProd(xit) + self.b) # Method B
            #yHat = sign(yit * (self.dotProd(xit) + self.b)) # Method C
            yHat = sign(yit * self.dotProd(xit)) # Method D

            ### if the value is actually good,
            ### and we predicted good, increment npr
            ### which is the # lines predicted as good that are actually good
            if yit == 1:
                self.testGood += 1
                if yHat > 0:
                    self.testnpr += 1
                    self.testdp += 1
                else:
                    self.testMistakes += 1
            if yit == -1:
                if yHat > 0:
                    self.testMistakes += 1

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
    def updateWeight(self, xarray, yStar):
        #startTime = time.time()

        ### update the weights in the weight vector that are in xit
        for element in xarray:
            self.w[element[0]] += self.eta * yStar * element[1]

        #endTime = time.time()
        #print 'updateWeight: ', endTime - startTime

    ### end updateWeight function
