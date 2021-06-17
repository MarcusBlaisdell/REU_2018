##########################################
### PEGASOS Test
### class
### Marcus Blaisdell
##########################################

import random
import time
import math
from functions import sign

### create a class for pegasos:
class pg_Class():
    ### global variables:
    trainData = []
    validationData = []
    testData = []
    T = 100 # default value of T
    k = 1 # default value of k
    lam = 1e-3 # lambda
    kmer = 11 # default k-mer size is 11
    #lamList = [5e-5, 5e-4, 5e-6, 1.29e-4, 3.07e-5,1.36e-4,1.67e-5,2e-4] # these are the lambda values from the Pegasos paper
    #lamList = [1e-5, 2e-5, 3e-5, 4e-5, 5e-5]
    lamList = [1e-3]
    #Tlist = [5, 10, 100, 1000, 5000, 10000]
    Tlist = [100]
    #klist = [1, 5, 10]
    klist = [1000, 2000, 5000, 10000]
    w = {} # weight is a dictionary
    eta = 0.0 # eta is calculated for each iteration
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
    nonZeroDict = {}

    def __init__(self):
        pass

    ### pegasos function:

    def pegasos(self, trainDataSet, testDataSet):
        # reset global mistake count so each iteration starts from zero
        self.trainMistakes = 0
        self.trainnpr = 0
        self.traindp = 0
        self.trainTotal = self.T*self.k
        xit = [] # x-sub-i-sub-t, the training vector for the current iterations
        yit = 0 # y-sub-i-sub-t, the label for the current training vector (yStar)
        yHat = 0

        ### repeat for T iterations:
        for t in range(self.T):
            ### Chose a subset from set trainData, uniformly at random
            ### the subset is contiguous, it is chosen at random from
            ### all subsets of S
            i = random.randint(0,len(trainDataSet) - self.k - 1)

            self.trainTotal += 1
            xit = trainDataSet[i][0]
            yit = trainDataSet[i][1]

            if yit == 1:
                self.trainGood += 1

            ### set eta = 1.0 / (lambda * t)
            ### (use t + 1 since the first iteration is t = 0)
            self.eta = 1.0 / (self.lam * (t + 1))
            ### both weight updates re-multiply eta by lambda,
            ### this is extra computation that isn't required,
            ### better to calculate the value once:
            tnew = 1.0 / float(t + 1)

            ### make the prediction:
            yHat = yit * self.dotProd(xit)

            ### update weight accordingly,
            ### if prediction is less than 1, update using trainData sample,
            ### if prediction is greater than or equal to 1, adjust weight
            ### by (1 - (1/t))
            if yHat < 1:
                ### since (1 - eta * lamda) is just factoring out
                ### the lambda from eta, replace by 'tnew'
                ### reduces the extra computation
                #self.w = (1 - self.eta * self.lam) * self.w + self.eta * yit * xit
                #self.w = (1 - tnew) * self.w + self.eta * yit * xit
                self.updateWeight(tnew, xit, yit)

                self.trainMistakes += 1
            else:
                #self.w = (1 - self.eta * self.lam) * self.w
                tScalar = 1 - tnew
                #self.w = (1 - tnew) * self.w
                ### only update non-zero indexes:
                nonZeroList = self.nonZeroDict.keys()
                for x in nonZeroList:
                    self.w[x] = tScalar * self.w[x]

                self.traindp += 1
                if yit == 1:
                    self.trainnpr += 1

        self.trainAccuracy = 100.0 - (100.0 * (self.trainMistakes / float(self.trainTotal)) )
        #print 'trainMistakes = ', self.trainMistakes, '   : train success = ', \
              #self.trainAccuracy, '%'
        if self.traindp > 0:
            self.trainPrecision = self.trainnpr / float(self.traindp)
        self.trainRecall = self.trainnpr / float(self.trainGood)
        #print 'trainPrecision: ', self.trainPrecision
        #print 'trainRecall: ', self.trainRecall
        if (self.trainPrecision + self.trainRecall) > 0:
            self.trainF1 = 2.0 * self.trainPrecision * self.trainRecall / (self.trainPrecision + self.trainRecall)
        #print 'F1: ', self.trainF1

    ### end Pegasos function

    ### new pegasosBatch function

    '''
    Initialize weight w to zero
    for t = 1, 2, ..., T:
        choose A_t as a proper subset of |m| where |A_t| = k; uniformly at random
        set A_t+ = {i in the set of A_t: yi<w_t, x_i> < 1}
        set eta_t = (1/(lambda * t))
        set w_t+1 = (1 - eta_t * lambda) * w_t + (eta_t / k) * sum i in the set of A_t+(y_i * x_i)
        w_t+1 = min {1, ((1/sqrt(lambda))norm(w_t+1))} * w_t+1

    '''

    def pegasosBatch (self, trainDataSet, testDataSet,  t):
        self.trainMistakes = 0
        self.trainnpr = 0
        self.traindp = 0

        # create an empty subset, A_t
        A_t = []
        # create an empty subset, A_t_plus
        A_t_plus = []

        # initialize w to zero:

        ### choose A_t as a proper subset of |m| where |A_t| = k; uniformly at random
        A_t = []
        for l in range(self.k):
            # select record uniformly at random
            i = random.randint(0, (len(trainDataSet) - 1) )
            A_t.append(trainDataSet[i])

        ### set A_t_plus = {i in the set of A_t: yi<w_t, x_i> < 1}
        A_t_plus = []
        for record in A_t:
            if ( (record[1] * (self.dotProd(record[0]) ) ) < 1 ):
                A_t_plus.append(record)
                self.trainMistakes += 1
            else:
                self.traindp += 1
                if record[1] == 1:
                    self.trainnpr += 1

        # set eta_t = (1/(lambda * t))
        self.eta = 1 / float(self.lam * (t + 1))

        # set w_t+1 = (1 - eta_t * lambda) * w_t + (eta_t / k) * sum i in the set of A_t+(y_i * x_i)

        ### sum i in the set of A_t+(y_i * x_i):

        ### runxit is a running dictionary of all unique indexes
        ### and the sums of their weights
        ### for each record in subset A_t_plus, sum the weights
        ### for each weight index:

        runxit = {}

        ### record[0] is the array of indexes and weights
        ### record[1] is the label of the vector, ({+1,-1})

        for record in A_t_plus:
            ### element[0] is the index of the weight
            ### element[1] is the weight of the vector at that index

            for element in record[0]:
                if runxit.get(element[0], '--') == '--':
                    runxit[element[0]] = element[1] * record[1]
                    self.w[element[0]] = 0
                else:
                    runxit[element[0]] += element[1] * record[1]

        ### convert runxit into an array of weight indexes and values:

        wArray = []

        wKeys = runxit.keys()
        for aKey in wKeys:
            wArray.append([aKey, runxit.get(aKey)])
            if self.nonZeroDict.get(aKey, '--') == '--':
                self.nonZeroDict[aKey] = 1

        ### set w_t+1 = (1 - eta_t * lambda) * w_t + (eta_t / k) * sum i in the set of A_t+(y_i * x_i)
        ### first, set w_t+1 = (1 - eta_t * lambda) * w_t:
        ### eta_t * lambda is just 1 / t so no need to use eta or lambda:
        tScalar = 1.0 - (1.0 / float(t + 1))
        nonZeroList = self.nonZeroDict.keys()
        for index in nonZeroList:
            self.w[index] = self.w[index] * tScalar
        #self.w = self.w * tScalar

        ### second, update weights of indexes from the sum times eta_t / k:
        eScalar = self.eta / float(self.k)
        for element in wArray:
            self.w[element[0]] += element[1] * eScalar

        ### Optional w_t+1 = min {1, ((1/sqrt(lambda))norm(w_t+1))} * w_t+1:
        #wNorm = LA.norm(self.w)
        wNorm = self.norm()
        if wNorm != 0:
            wOpt = 1.0 / float(math.sqrt(self.lam) * wNorm)
        else:
            wOpt = 0

        if wOpt < 1:
            nonZeroList = self.nonZeroDict.keys()
            for index in nonZeroList:
                self.w[index] = self.w[index] * wOpt
            #self.w = self.w * wOpt

        self.trainAccuracy = 100.0 - (100.0 * (self.trainMistakes / float(((t + 1) * self.k))))
        #print 'trainAccuracy: ', self.trainAccuracy

    ### end new pegasosBatch function

    ### use new weight to test accuracy on test dataList

    def testWeight (self, trainDataSet, testDataSet):
        # reset mistakes count so each iteration starts at 0
        self.testMistakes = 0
        self.testnpr = 0
        self.testdp = 0
        self.testTotal = len(testDataSet)
        xit = [] # x-sub-i-sub-t, the training vector for the current iterations
        yit = 0 # y-sub-i-sub-t, the label for the current training vector (yStar)
        yHat = 0

        aVal = 0
        bVal = 0
        cVal = 0
        dVal = 0
        ### evaluate all test samples:
        for i in range(len(testDataSet)):
            xit = testDataSet[i][0]
            yit = testDataSet[i][1]

            ### make the prediction:
            #yHat = yit * self.dotProdTest(xit)
            yHat = sign(yit * self.dotProd(xit))
            #print 'self.dotProdTest(xit): ', self.dotProdTest(xit), ' :'
            #print 'yit * self.dotProdTest(xit)', yit * self.dotProdTest(xit)
            #print 'yHat: ', yHat, 'yit: ', yit


            ### if the value is actually good,
            ### and we predicted good, increment npr
            ### which is the # lines predicted as good that are actually good
            if yit == 1:
                if yHat >= 1:
                    self.testnpr += 1
                else:
                    self.testMistakes += 1
            else: #if yit == -1:
                if yHat >= 1:
                    self.testMistakes += 1

            ### evaluate the prediction
            ### if predicted good, increment predicted good count, dp
            if yHat >= 1:
                self.testdp += 1

        self.testAccuracy = 100.0 - (100.0 * (self.testMistakes / float(self.testTotal)) )
        #print 'testMistakes = ', self.testMistakes, '   : test success = ', \
              #self.testAccuracy, '%'

        if self.testdp > 0:
            self.testPrecision = self.testnpr / float(self.testdp)
        self.testRecall = self.testnpr / float(self.testGood)
        #print 'testPrecision: ', self.testPrecision
        #print 'testRecall: ', self.testRecall
        if (self.testPrecision + self.testRecall) > 0:
            self.testF1 = 2.0 * self.testPrecision * self.testRecall / (self.testPrecision + self.testRecall)
    ### end testWeight function

    ### get the count of the total number of train / test samples
    ### that are actually good:
    def countGood (self, trainDataSet, testDataSet):
        for i in range(len(trainDataSet)):
            if trainDataSet[i][1] == 1:
                self.trainGood += 1

        for i in range(len(self.validationData)):
            if self.validationData[i][1] == 1:
                self.validationGood += 1

        for i in range(len(testDataSet)):
            if testDataSet[i][1] == 1:
                self.testGood += 1

    ### end countGood function


    ### function dotProd(), dot product by index
    def dotProd(self, xArray):
        #startTime = time.time()
        runSum = 0.0

        ### xArray format is [[w[index],value],...]
        ### xArray is a sparse vector, anything that is not in it will be zero,
        ### so only sum the products of indexes from w that exist in xArray
        for element in xArray:
            if self.w.get(element[0], '--') != '--':
                runSum += self.w[element[0]] * element[1]

        #endTime = time.time()
        #print 'dotProd2: ', endTime - startTime

        return runSum

    ### end dotProd() function

    ### function dotProdTest(), dot product by index
    def dotProdTest(self, xArray):
        #startTime = time.time()
        result = 0.0

        ### xArray format is [[w[index],value],...]
        ### xArray is a sparse vector, anything that is not in it will be zero,
        ### so only sum the products of indexes from w that exist in xArray
        for element in xArray:
            if self.w.get(element[0], '--') != '--':
                result += self.w[element[0]] * element[1]

        #endTime = time.time()
        #print 'dotProd2: ', endTime - startTime

        #print 'result: ', result
        return result

    ### end dotProdTest() function

    ### updateWeight function:

    def updateWeight(self, tnew, xarray, yStar):
        #startTime = time.time()
        tScalar = 1 - tnew
        wOpt = 0.0
        #wNorm = LA.norm(self.w)
        wNorm = self.norm()
        wOpt = wNorm / math.sqrt(self.lam)
        if wOpt > 1 or math.isnan(wOpt) or math.isinf(wOpt):
            wOpt = 1

        ### This only processes non-zero features:
        keyList = self.nonZeroDict.keys()
        for key in keyList:
            self.w[key] = self.w[key] * tScalar
        #endTime = time.time()
        #print 'iterationWeight: ', endTime - startTime

        #startTime = time.time()
        for element in xarray:
            self.w[element[0]] += self.eta * yStar * element[1]
            ### if a given index of the weight vector is non-zero,
            ### add it to the dictionary if it isn't already in it
            if self.w[element[0]] != 0:
                if self.nonZeroDict.get(element[0], '-') == '-':
                    self.nonZeroDict[element[0]] = 1
                    self.w[element[0]] = 0
            ### if the weight value at that index has become a zero,
            ### delete it from the non-zero list
            else:
                if self.nonZeroDict.get(element[0], '-') != '-':
                    del self.nonZeroDict[element[0]]

        ### optional scalar:
        for key in keyList:
            self.w[key] = self.w[key] * wOpt
        '''
            ### This is slower than above by a notable amount
            ### average runtime 22.57s vs. 19.99s = 2.58s longer
            if self.nonZeroDict.get(element[0], '-') == '-':
                self.nonZeroDict[element[0]] = 1
        '''

        #endTime = time.time()
        #print 'updateWeight: ', endTime - startTime

    ### end updateWeight function

    ### updateWeightBatch function:

    def updateWeightBatch(self, tnew, xarray, yStar):
        tScalar = 1 - tnew
        #startTime = time.time()

        ### This only processes non-zero features:

        keyList = self.nonZeroDict.keys()
        for key in keyList:
            self.w[key] = self.w[key] * tScalar
        #endTime = time.time()
        #print 'iterationWeight: ', endTime - startTime

        #startTime = time.time()
        for element in xarray:
            self.w[element[0]] += self.eta / self.k * yStar * element[1]
            ### if a given index of the weight vector is non-zero,
            ### add it to the dictionary if it isn't already in it
            if self.w[element[0]] != 0:
                if self.nonZeroDict.get(element[0], '-') == '-':
                    self.nonZeroDict[element[0]] = 1
                    self.w[element[0]] = 0
            ### if the weight value at that index has become a zero,
            ### delete it from the non-zero list
            else:
                if self.nonZeroDict.get(element[0], '-') != '-':
                    del self.nonZeroDict[element[0]]

        '''
        ### optional scalar:

        wOpt = 0.0
        wNorm = LA.norm(self.w) # norm of current weight vector
        wOpt = 1.0 / (wNorm * math.sqrt(self.lam) )
        if wOpt > 1 or math.isnan(wOpt) or math.isinf(wOpt):
            wOpt = 1

        for key in keyList:
            self.w[key] = self.w[key] * wOpt

        ### end optional scalar
        '''

        #endTime = time.time()
        #print 'updateWeight: ', endTime - startTime

    ### end updateWeightBatch function

    ### norm function:
    def norm(self):
        runSum = 0.0
        wIndex = self.w.keys()
        for i in wIndex:
            runSum += self.w.get(i)**2
        return math.sqrt(runSum)
