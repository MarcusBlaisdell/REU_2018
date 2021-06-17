##########################################
### Perceptron Test
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
    T = 10 # default value of T
    k = 1 # default value of k
    kmer = 11 # default k-mer size is 11
    #Tlist = [5, 10, 100, 1000, 5000, 10000]
    Tlist = [20]
    #klist = [1, 5, 10]
    klist = [1]
    w = {} # weight is a dictionary
    eta = 1.0 # eta is learning rate
    lam = 1e-15
    l1Norm = 0.0
    underSample = 0
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
        #print 'perceptron called'
        xit = [] # x-sub-i-sub-t, the training vector for the current iterations
        yit = 0 # y-sub-i-sub-t, the label for the current training vector (yStar)
        yHat = 0

        for l in range (len(trainDataSet)):
            self.trainTotal += 1
            xit = trainDataSet[l][0]
            yit = trainDataSet[l][1]
            if yit == 1:
                self.trainGood += 1 ### increment count of actual good

            ### make the prediction: yHat = y*(<w,x>)
            #yHat = yit * self.dotProd(xit) # Method U
            #yHat = yit * (self.dotProd(xit) + self.b) # Method V
            #yHat = sign(yit * (self.dotProd(xit) + self.b)) # Method W
            yHat = sign(yit * self.dotProd(xit)) # Method X

            ### l1 regularization
            #yHat = yit * self.dotProd(xit) + self.lam * self.l1Norm
            '''
            ### hinge loss:
            yHat = yit * self.dotProd (xit)
            if yHat >= 0:
                yHat = 0
            '''

            #print 'initial eval: ', yit, ' : ', yHat

            ### update weight accordingly,
            ### if predicted value and actual value don't match,
            ## update weight,
            ### if they do match no update required

            if yHat > 0:
                self.traindp += 1 # +1 prediction count
                if yit == 1:
                    self.trainnpr += 1 # correct +1 prediction
            else: ### yHat <= 0:
                self.trainMistakes += 1
                self.updateWeight(xit, yit, yHat)
                if yit == 1:
                    self.traindp += 1 # +1 prediction count
            '''
            if yHat <= 0:
                self.trainMistakes += 1
                self.updateWeight(xit, yit)
                if yit == 1:
                    self.traindp += 1 # +1 prediction count
            else:
                self.traindp += 1 # +1 prediction count
                if yit == 1:
                    self.trainnpr += 1 # correct +1 prediction
            '''

            '''
            if yit == 1:
                self.trainGood += 1
                if yHat > 0: #prediction is +1:
                    self.trainnpr += 1 # correct +1 prediction
                    self.traindp += 1 # +1 prediction count
                else: # prediction is -1:
                    #print yit, ' : ', yHat
                    self.updateWeight(xit, yit)
                    self.trainMistakes += 1
            else: # yit == -1
                if yHat > 0: # if prediction is +1:
                    #print yit, ' : ', yHat
                    self.traindp += 1 # +1 prediction count:
                    self.trainMistakes += 1
                    self.updateWeight(xit, yit)
            '''

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
            if yit == 1:
                self.testGood += 1 ### increment count of actual good

            ### make the prediction:
            #yHat = yit * self.dotProd(xit) # Method A
            #yHat = yit * (self.dotProd(xit) + self.b) # Method B
            #yHat = sign(yit * (self.dotProd(xit) + self.b)) # Method C
            yHat = sign(yit * self.dotProd(xit)) # Method D

            ### if the value is actually good,
            ### and we predicted good, increment npr
            ### which is the # lines predicted as good that are actually good
            if yHat > 0:
                self.testdp += 1 # +1 prediction count
                if yit == 1:
                    self.testnpr += 1 # correct +1 prediction
            else: ### yHat <= 0:
                self.testMistakes += 1
                if yit == 1:
                    self.testdp += 1 # +1 prediction count
            '''
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
                    self.testdp += 1
            '''

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
            '''
            ### the l1 norm can be maintained dynamically here
            ### so it doesn't have to be re-computed at each evaluation

            curW = self.w[element[0]]
            newW = curW + self.eta * yStar * element[1]
            self.l1Norm = self.l1Norm - abs(curW) + abs(newW)
            '''

            ### compensate for imbalanced data:
            ### if it should be good, use a larger learning rate
            ### if it should be bad, use a smaller one:

            if yStar == 1 and yHat == -1:
                self.eta = 2.0
            else:
                self.eta = 0.2
            '''
            ### under sample irrelvant records by 1/4
            if yStar == 1:
                if self.underSample == 15:
                    self.w[element[0]] += self.eta * yStar * element[1]
                    self.underSample = 0
                else:
                    self.underSample += 1
            else:
                self.w[element[0]] += self.eta * yStar * element[1]
            '''

            self.w[element[0]] += self.eta * yStar * element[1]
            #print 'weight: ', self.w[element[0]], '   yStar: ', yStar, '   element[1]', element[1]

        #endTime = time.time()
        #print 'updateWeight: ', endTime - startTime

    ### end updateWeight function

    ### theNorm function
    ### returns the l1 norm of the weight vector

    def theNorm (self):
        runningSum = 0.0

        indexList = self.w.keys()

        for element in indexList:
            runningSum += self.w.get(element)

        return runningSum

    ### end theNorm function
