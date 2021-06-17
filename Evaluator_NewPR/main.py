from sklearn.metrics import precision_score
import numpy as np
from pr_Class import *
from functions import *
import logging
import sys
import time

#geneList = ['sonarTrainBalanced']
#geneTest = 'sonarTest'
#geneList = ['smallTrain']
#geneTest = 'smallTest'
geneList = ['2016C', 'CH611', 'Co6114']
#geneTest = '2016C'
geneTest = '0157'
kmerList = [11]
kmer = 11
maxReadSize = 10000

thePath = "/Users/MarcusBlaisdell/Documents/LinuxShare/tenK/"
resultsStringPre = "results_perceptron_P-at-R-98_" + geneTest + "_"
resultsStringPost = ".csv"

### create empty log file:
logging.basicConfig(filename='run.log.perceptron_P-at-R-98_',level=logging.DEBUG)
logging.info("Log File created")

myPRClass = pr_Class()

myPRClass.kmer = kmer

resultsFile = resultsStringPre + str(kmer) + resultsStringPost
print 'resultsFile: ', resultsFile
logging.info ('results File: ' + resultsFile)

print 'set: ', str(kmer)
logging.info ('set: ' + str(str(kmer)))

trainDataSet = []
testDataSet = []

print 'opening results file:'
logging.info ('opening results file:' + '\n')
outFile = open(resultsFile, "w")

### Begin Perceptron:

outFile.write ('Perceptron'  + '\n' + '\n')
outFile.write ('Gene left out: ' + geneTest + '\n')
outFile.write('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
              'Accuracy' + ',' + 'T' + ',' + 'K' + ',' \
              + 'Mistakes' + ',' + 'Iteration' + '\n')
logging.info('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
              'Accuracy' + ',' + 'T' + ',' + 'K' + ',' \
              + 'Mistakes' + ',' + 'Iteration' + '\n')

print 'running perceptron . . . '
logging.info ('running perceptron . . . ' + '\n')

startTime = time.time()

for t in range(myPRClass.T):
    myPRClass.trainnpr = 0
    myPRClass.traindp = 0
    myPRClass.trainMistakes = 0
    myPRClass.trainTotal = 0
    myPRClass.trainGood = 0

    myPRClass.testnpr = 0
    myPRClass.testdp = 0
    myPRClass.testMistakes = 0
    myPRClass.testTotal = 0
    myPRClass.testGood = 0

    myPRClass.scores = []
    myPRClass.labels = []

    for gene in geneList:
        trainTime1 = time.time()
        readGene (thePath, gene, kmer, maxReadSize, "perceptron", myPRClass, t)
        trainTime = time.time() - trainTime1
        outFile.write(str(trainTime) + '\n')
        logging.info (str(trainTime) + '\n')

        ### report results of training

        myPRClass.trainAccuracy = 100 - (100 * (myPRClass.trainMistakes / float(myPRClass.trainTotal)) )

        if myPRClass.traindp > 0:
            myPRClass.trainPrecision = myPRClass.trainnpr / float(myPRClass.traindp)
        if myPRClass.trainGood > 0:
            myPRClass.trainRecall = myPRClass.trainnpr / float(myPRClass.trainGood)

        if (myPRClass.trainPrecision + myPRClass.trainRecall) > 0:
            myPRClass.trainF1 = 2 * myPRClass.trainPrecision * myPRClass.trainRecall / (myPRClass.trainPrecision + myPRClass.trainRecall)

        outFile.write('train' + ',' + str(myPRClass.trainPrecision) + ',' + str(myPRClass.trainRecall) \
                      + ',' + str(myPRClass.trainF1) + ',' + str(myPRClass.trainAccuracy)\
                      + ',' + str(myPRClass.T) + ',' + str(myPRClass.k) + ',' \
                      + str(myPRClass.trainMistakes) + ',' + str(t + 1) + '\n')
        logging.info ('train' + ',' + str(myPRClass.trainPrecision) + ',' + str(myPRClass.trainRecall) \
                      + ',' + str(myPRClass.trainF1) + ',' + str(myPRClass.trainAccuracy)\
                      + ',' + str(myPRClass.T) + ',' + str(myPRClass.k) + ',' \
                      + str(myPRClass.trainMistakes) + ',' + str(t + 1) + '\n')

        logging.info('training Complete')
        #####

    ### Test

    testTime1 = time.time()
    readGene (thePath, geneTest, kmer, maxReadSize, "testWeight", myPRClass, t)
    testTime = time.time() - testTime1
    outFile.write (str(testTime) + '\n')
    logging.info (str(testTime) + '\n')

    scores = np.array(myPRClass.scores)
    labels = np.array(myPRClass.labels)
    positive_scores = scores[labels == 1]
    print 'len(positive_scores): ', len(positive_scores)
    threshold = np.percentile (positive_scores, myPRClass.target_recall * 100)
    predicted = scores >= threshold
    print 'len(predicted): ', len(predicted)
    print 'predicted: ', predicted
    print 'labels: ', labels
    pScore = precision_score (labels, predicted)
    print 'precision score at ', myPRClass.target_recall, ': ', pScore
    outFile.write ('precision score at ' + str(myPRClass.target_recall) + ' recall: ' + str(pScore) + '\n')
    logging.info ('precision score at ' + str(myPRClass.target_recall) + ' recall: ' + str(pScore) + '\n')

    ### report results of testing

    myPRClass.testAccuracy = 100 - (100 * (myPRClass.testMistakes / float(myPRClass.testTotal)) )
    #print 'testMistakes = ', myPRClass.testMistakes, '   : test success = ', \
          #myPRClass.testAccuracy, '%'
    if myPRClass.testdp > 0:
        myPRClass.testPrecision = myPRClass.testnpr / float(myPRClass.testdp)
    if myPRClass.testGood > 0:
        myPRClass.testRecall = myPRClass.testnpr / float(myPRClass.testGood)
    #print 'testPrecision: ', myPRClass.testPrecision
    #print 'testRecall: ', myPRClass.testRecall
    if (myPRClass.testPrecision + myPRClass.testRecall) > 0:
        myPRClass.testF1 = 2 * myPRClass.testPrecision * myPRClass.testRecall / (myPRClass.testPrecision + myPRClass.testRecall)
    #print 'F1: ', myPRClass.testF1
    outFile.write('test' + ',' + str(myPRClass.testPrecision) + ',' + str(myPRClass.testRecall) \
                  + ',' + str(myPRClass.testF1) + ',' + str(myPRClass.testAccuracy)\
                  + ',' + str(myPRClass.T) + ',' + str(myPRClass.k) + ',' \
                  + str(myPRClass.testMistakes) + ',' + str(t + 1) + '\n')
    logging.info ('test' + ',' + str(myPRClass.testPrecision) + ',' + str(myPRClass.testRecall) \
                  + ',' + str(myPRClass.testF1) + ',' + str(myPRClass.testAccuracy)\
                  + ',' + str(myPRClass.T) + ',' + str(myPRClass.k) + ',' \
                  + str(myPRClass.testMistakes) + ',' + str(t + 1) + '\n')
    #####
