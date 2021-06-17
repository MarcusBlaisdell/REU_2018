##########################################
### Perceptron Test
### class
### Marcus Blaisdell
##########################################

import random
import time
import numpy as np

class p_Class():
    ### global variables:
    trainData = []
    validationData = []
    testData = []
    T = 20 # default value of T
    k = 1 # default value of k
    kmer = 11 # default k-mer size is 11
    #Tlist = [5, 10, 100, 1000, 5000, 10000]
    Tlist = [20]
    #klist = [1, 5, 10]
    klist = [1]
    w = np.array([]) # weight is an numpy array
    b = 0.0 # bias term
    eta = 1.0 # eta is learning rate
    lam = 1e-1 # regularization parameter
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

    ### loadDatafunction

    def loadData(self, fileName, dataName):
        word = ''
        value = ''
        line = []
        label = 0
        readFile = open(fileName, "r")

        stringOne = readFile.readline()

        if stringOne[0] == 'f':
            ### parse out header from first line:
            i = 4
            flag = 0
            while flag != 1:
                if stringOne[i] == 'a' and stringOne[i-1] == 'p' and stringOne[i-2] == 'p' and stringOne[i-3] == 'o':
                    flag = 1
                else:
                    i += 1

            flag = 0
            while flag == 0:
                try:
                    if int (stringOne[i]) > 0:
                        flag = 1
                    else:
                        i+=1
                except:
                    i+=1

            string = stringOne[i:]
            ### end parse Header
        else:
            string = stringOne

        ### now, parse data into a list
        while string:
            i = 0

            if string[len(string) - 3] == '-':
                label = -1
            else:
                label = 1

            while i < len(string):
                if string[i] == ':':
                    i += 1
                    while string[i] != ' ':
                        value += string[i]
                        i += 1
                    tuple = [int(word), float(value)]
                    line.append(tuple)
                    word = ''
                    value = ''

                else:
                    word += string[i]
                i += 1

            lineTuple = [line, label]
            #self.trainData.append(lineTuple)
            dataName.append(lineTuple)
            line = []
            word = ''
            string = readFile.readline()

        readFile.close()

    ### end loadData functions

    ### perceptron function:

    def perceptron(self, outFile):
        self.trainTotal = len(self.trainData)
        xit = np.array([]) # x-sub-i-sub-t, the training vector for the current iterations
        yit = 0 # y-sub-i-sub-t, the label for the current training vector (yStar)
        yHat = 0

        ### initialize weight vector:
        #self.w = np.zeros(len(self.trainData[0][0]), dtype=float)
        ### weight vector size is {[4^(len(k-mer)) - len(k-mer)] + 1}

        kmerSize = ((4**self.kmer) + 1)
        self.w = np.zeros(kmerSize, dtype=float)
        #self.w = np.zeros(10000000, dtype=float)

        ### repeat for T iterations:
        for t in range(self.T):
            # reset global mistake count so each iteration starts from zero
            self.trainMistakes = 0
            self.trainnpr = 0
            self.traindp = 0
            #print 'Iteration: ', t
            for l in range (len(self.trainData)):
                xit = np.array(self.trainData[l][0])
                yit = self.trainData[l][1]

                ### make the prediction: yHat = y*(<w,x>)
                yHat = yit * self.dotProd(xit)
                #print 'initial eval: ', yit, ' : ', yHat

                ### update weight accordingly,
                ### if predicted value and actual value don't match,
                ## update weight,
                ### if they do match no update required
                if yit == 1:
                    if yHat > 0:
                        self.trainnpr += 1
                        self.traindp += 1
                    else:
                        #print yit, ' : ', yHat
                        self.updateWeight(xit, yit)
                        self.trainMistakes += 1
                else:
                    if yHat <= 0:
                        #print yit, ' : ', yHat
                        self.traindp += 1
                        self.trainMistakes += 1
                        self.updateWeight(xit, yit)


                '''
                if (yHat != yit):
                    self.updateWeight(xit, yit)
                    ### update incorrect prediction count
                    self.trainMistakes += 1
                else:
                    ### if prediction was correct, update correct prediction counts
                    self.traindp += 1
                    self.trainnpr += 1
                '''

                ### end train loop, trains on each sample in trainData

            ### report results of training

            self.trainAccuracy = 100 - (100 * (self.trainMistakes / float(self.trainTotal)) )
            #print 'trainMistakes = ', self.trainMistakes, '   : train success = ', \
                  #self.trainAccuracy, '%'
            if self.traindp > 0:
                self.trainPrecision = self.trainnpr / float(self.traindp)
            self.trainRecall = self.trainnpr / float(self.trainGood)
            #print 'trainPrecision: ', self.trainPrecision
            #print 'trainRecall: ', self.trainRecall
            if (self.trainPrecision + self.trainRecall) > 0:
                self.trainF1 = 2 * self.trainPrecision * self.trainRecall / (self.trainPrecision + self.trainRecall)
            #print 'F1: ', self.trainF1
            outFile.write('train' + ',' + str(self.trainPrecision) + ',' + str(self.trainRecall) \
                          + ',' + str(self.trainF1) + ',' + str(self.trainAccuracy)\
                          + ',' + str(self.T) + ',' + str(self.k) + ',' \
                          + str(self.trainMistakes) + ',' + str(t + 1) + '\n')

            #self.validationWeight(outFile, t)
            self.testWeight(outFile, t)

        ### end iteration loop

    ### end perceptron function

    ### perceptronBias function:

    def perceptronBias(self, outFile):
        self.trainTotal = len(self.trainData)
        xit = np.array([]) # x-sub-i-sub-t, the training vector for the current iterations
        yit = 0 # y-sub-i-sub-t, the label for the current training vector (yStar)
        yHat = 0

        ### initialize weight vector:
        #self.w = np.zeros(len(self.trainData[0][0]), dtype=float)
        ### weight vector size is {[4^(len(k-mer)) - len(k-mer)] + 1}

        kmerSize = ((4**self.kmer) + 1)
        self.w = np.zeros(kmerSize, dtype=float)
        #self.w = np.zeros(10000000, dtype=float)

        ### repeat for T iterations:
        for t in range(self.T):
            # reset global mistake count so each iteration starts from zero
            self.trainMistakes = 0
            self.trainnpr = 0
            self.traindp = 0
            #print 'Iteration: ', t
            for l in range (len(self.trainData)):
                xit = np.array(self.trainData[l][0])
                yit = self.trainData[l][1]

                ### make the prediction: yHat = y*(<w,x>)
                yHat = yit * (self.dotProd(xit) + self.b)
                #print 'initial eval: ', yit, ' : ', yHat

                ### update weight accordingly,
                ### if predicted value and actual value don't match,
                ## update weight,
                ### if they do match no update required
                if yit == 1:
                    if yHat > 0:
                        self.trainnpr += 1
                        self.traindp += 1
                    else:
                        #print yit, ' : ', yHat
                        self.updateWeight(xit, yit)
                        self.trainMistakes += 1
                else:
                    if yHat <= 0:
                        #print yit, ' : ', yHat
                        self.traindp += 1
                        self.trainMistakes += 1
                        self.updateWeight(xit, yit)


                '''
                if (yHat != yit):
                    self.updateWeight(xit, yit)
                    ### update incorrect prediction count
                    self.trainMistakes += 1
                else:
                    ### if prediction was correct, update correct prediction counts
                    self.traindp += 1
                    self.trainnpr += 1
                '''

                ### end train loop, trains on each sample in trainData

            ### report results of training

            self.trainAccuracy = 100 - (100 * (self.trainMistakes / float(self.trainTotal)) )
            #print 'trainMistakes = ', self.trainMistakes, '   : train success = ', \
                  #self.trainAccuracy, '%'
            if self.traindp > 0:
                self.trainPrecision = self.trainnpr / float(self.traindp)
            self.trainRecall = self.trainnpr / float(self.trainGood)
            #print 'trainPrecision: ', self.trainPrecision
            #print 'trainRecall: ', self.trainRecall
            if (self.trainPrecision + self.trainRecall) > 0:
                self.trainF1 = 2 * self.trainPrecision * self.trainRecall / (self.trainPrecision + self.trainRecall)
            #print 'F1: ', self.trainF1
            outFile.write('train' + ',' + str(self.trainPrecision) + ',' + str(self.trainRecall) \
                          + ',' + str(self.trainF1) + ',' + str(self.trainAccuracy)\
                          + ',' + str(self.T) + ',' + str(self.k) + ',' \
                          + str(self.trainMistakes) + ',' + str(t + 1) + '\n')

            #self.validationWeight(outFile, t)
            self.testWeight(outFile, t)

        ### end iteration loop

    ### end perceptronBias function

    ### use new weight to test accuracy on validation dataList

    def validationWeight (self, outFile, t):
        # reset mistakes count so each iteration starts at 0
        self.validationMistakes = 0
        self.validationnpr = 0
        self.validationdp = 0
        self.validationTotal = len(self.validationData)
        xit = np.array([]) # x-sub-i-sub-t, the training vector for the current iterations
        yit = 0 # y-sub-i-sub-t, the label for the current training vector (yStar)
        yHat = 0

        ### evaluate all test samples:
        for i in range(len(self.validationData)):
            xit = np.array(self.validationData[i][0])
            yit = self.validationData[i][1]

            ### make the prediction:
            #yHat = yit * np.dot(self.w, xit)
            yHat = yit * self.dotProd(xit)

            ###
            ### if the value is actually good,
            ### and we predicted good, increment npr
            ### which is the # lines predicted as good that are actually good
            if yit == 1:
                if yHat > 0:
                    self.validationnpr += 1
                    self.validationdp += 1
                else:
                    self.validationMistakes += 1
            if yit == -1:
                if yHat <= 0:
                    self.validationMistakes += 1
                    self.validationdp += 1

            '''
            ### evaluate the prediction
            ### if predicted good, increment predicted good count, dp
            if yHat == 1:
                self.validationdp += 1
            '''

        self.validationAccuracy = 100 - (100 * (self.validationMistakes / float(self.validationTotal)) )
        #print 'validationMistakes = ', self.validationMistakes, '   : validation success = ', \
              #self.validationAccuracy, '%'
        if self.validationdp > 0:
            self.validationPrecision = self.validationnpr / float(self.validationdp)
        self.validationRecall = self.validationnpr / float(self.validationGood)
        #print 'validationPrecision: ', self.validationPrecision
        #print 'validationRecall: ', self.validationRecall
        if (self.validationPrecision + self.validationRecall) > 0:
            self.validationF1 = 2 * self.validationPrecision * self.validationRecall / (self.validationPrecision + self.validationRecall)

        outFile.write('validation' + ',' + str(self.validationPrecision) + ',' + str(self.validationRecall) \
                      + ',' + str(self.validationF1) + ',' + str(self.validationAccuracy)\
                      + ',' + str(self.T) + ',' + str(self.k) + ',' \
                      + str(self.validationMistakes) + ',' + str(t + 1) + '\n')

    ### end validationWeight function

    ### use new weight to test accuracy on test dataList

    def testWeight (self, outFile, t):
        # reset mistakes count so each iteration starts at 0
        self.testMistakes = 0
        self.testnpr = 0
        self.testdp = 0
        self.testTotal = len(self.testData)
        xit = np.array([]) # x-sub-i-sub-t, the training vector for the current iterations
        yit = 0 # y-sub-i-sub-t, the label for the current training vector (yStar)
        yHat = 0

        ### evaluate all test samples:
        for i in range(len(self.testData)):
            xit = np.array(self.testData[i][0])
            yit = self.testData[i][1]

            ### make the prediction:
            #yHat = yit * np.dot(self.w, xit)
            yHat = yit * self.dotProd(xit)

            ### if the value is actually good,
            ### and we predicted good, increment npr
            ### which is the # lines predicted as good that are actually good
            if yit == 1:
                if yHat > 0:
                    self.testnpr += 1
                    self.testdp += 1
                else:
                    self.testMistakes += 1
            if yit == -1:
                if yHat <= 0:
                    self.testMistakes += 1
                    self.testdp += 1

            '''
            ### evaluate the prediction
            ### if predicted good, increment predicted good count, dp
            if yHat == 1:
                self.testdp += 1
            '''

        self.testAccuracy = 100 - (100 * (self.testMistakes / float(self.testTotal)) )
        #print 'testMistakes = ', self.testMistakes, '   : test success = ', \
              #self.testAccuracy, '%'

        if self.testdp > 0:
            self.testPrecision = self.testnpr / float(self.testdp)
        self.testRecall = self.testnpr / float(self.testGood)
        #print 'testPrecision: ', self.testPrecision
        #print 'testRecall: ', self.testRecall
        if (self.testPrecision + self.testRecall) > 0:
            self.testF1 = 2 * self.testPrecision * self.testRecall / (self.testPrecision + self.testRecall)

        outFile.write('test' + ',' + str(self.testPrecision) + ',' + str(self.testRecall) \
                      + ',' + str(self.testF1) + ',' + str(self.testAccuracy)\
                      + ',' + str(self.T) + ',' + str(self.k) + ',' \
                      + str(self.testMistakes) + ',' + str(t + 1) + '\n')

    ### end testWeight function

    ### get the count of the total number of train / test samples
    ### that are actually good:
    def countGood (self):
        for i in range(len(self.trainData)):
            if self.trainData[i][1] == 1:
                self.trainGood += 1

        for i in range(len(self.validationData)):
            if self.validationData[i][1] == 1:
                self.validationGood += 1

        for i in range(len(self.testData)):
            if self.testData[i][1] == 1:
                self.testGood += 1

    ### end countGood function

    '''
    ### function dotProd(), dot product by index
    def dotProd(self, xArray):
        #startTime = time.time()
        result = 0.0

        for element in xArray:
            index = element[0]
            value = element[1]

            result += self.w[index] * value

        #endTime = time.time()
        #print 'dotProd: ', endTime - startTime

        return result

    ### end dotProd() function
    '''
    ### function dotProd(), dot product by index
    def dotProd(self, xArray):
        #startTime = time.time()
        result = 0.0

        for element in xArray:
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
            self.w[element[0]] += self.eta * self.lam * yStar * element[1]

        ### update the bias term:
        #self.b = self.b +  100 * yStar

        #endTime = time.time()
        #print 'updateWeight: ', endTime - startTime

    ### end updateWeight function

    ### sign function
    def sign(self, number):
        if number > 0:
            return 1
        else:
            return -1

    ### end sign function
