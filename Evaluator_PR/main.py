##########################################
### PEGASOS versus Perceptron Test
### for genomic data
### Marcus Blaisdell
### Parallel Processing
###
##########################################

import multiprocessing as mp
from pg_Class import *
from pr_Class import *
from functions import *
import logging
import sys
import time

### create array of the values to label the training, testing, results files:

#geneList = ['0157', '2016C', 'CH611', 'Co6114', 'ED1a', 'EDL933-1', 'FAP1', '_isolate102', 'RS76', 'UMN026']
#geneList = ['2016C', 'CH611', 'Co6114']
#geneList = ['0157', 'Co6114', '2016C']
geneList = ['sonarTrainBalanced']
#geneList = ['2016C']
#geneTest = '0157'
#geneTest = 'CH611'
geneTest = 'sonarTest'
#kmerList = [11, 13, 15, 17]
kmerList = [11]
#subList = ['a', 'b', 'c', 'd']
subList = ['a']
#biasVariableList = [1, 10, 50, 100]
biasVariableList = [0]

maxReadSize = 10000

### run the comparison for each training/testing set:

#thePath = "/home/marcus/Data/Genomic-Data/Testing/"
#thePath = "/data/doppa/users/mblaisdell/genome-assembly-datasets/"
thePath = "/Users/MarcusBlaisdell/Documents/LinuxShare/tenK/"
#thePath = ""
#resultsStringPre = "/home/marcus/Documents/LinuxShare/Results/resultszU-A-b-1-100_"
resultsStringPre = "results_perceptron_corrected_sonarb_" + geneTest + "_"
resultsStringPost = ".csv"

### create empty log file:
logging.basicConfig(filename='run.log.perceptron3',level=logging.DEBUG)
logging.info("Log File created")

def individualProcess(kmer, letter):
    myPGClass = pg_Class()
    myPRClass = pr_Class()

    ### The datasets are the same for each function, create them once to save memmory:
    trainDataSet = []
    testDataSet = []

    myPGClass.kmer = kmer
    myPRClass.kmer = kmer

    resultsFile = resultsStringPre + str(kmer) + resultsStringPost
    print 'resultsFile: ', resultsFile
    logging.info ('results File: ' + resultsFile)

    print 'set: ', str(kmer)
    logging.info ('set: ' + str(str(kmer)))

    print 'loading data:'
    logging.info('loading data')

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
                  + 'Bias' + ','\
                  + 'Mistakes' + ',' + 'Iteration' + '\n')
    logging.info('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
                  'Accuracy' + ',' + 'T' + ',' + 'K' + ',' \
                  + 'Bias' + ','\
                  + 'Mistakes' + ',' + 'Iteration' + '\n')

    print 'running perceptron . . . '
    logging.info ('running perceptron . . . ' + '\n')

    startTime = time.time()
    for biasVariable in biasVariableList:
        myPRClass.b = biasVariable
        #####
        for t in range(myPRClass.T):
            #myPRClass.eta = 1 / (t + 1)
            #print 'iteration: ', t
            #####
            myPRClass.trainnpr = 0
            myPRClass.traindp = 0
            myPRClass.trainMistakes = 0
            myPRClass.trainTotal = 0
            myPRClass.trainGood = 0

            for gene in geneList:
                trainTime1 = time.time()
                readGene (thePath, gene, kmer, maxReadSize, "perceptron", myPRClass, t)
                trainTime = time.time() - trainTime1
                outFile.write ('trainTime: ' + str(trainTime))
                logging.info ('trainTime: ' + str(trainTime))

                ### report results of training

                myPRClass.trainAccuracy = 100 - (100 * (myPRClass.trainMistakes / float(myPRClass.trainTotal)) )
                #print 'trainMistakes = ', myPRClass.trainMistakes, '   : train success = ', \
                      #myPRClass.trainAccuracy, '%'
                if myPRClass.traindp > 0:
                    myPRClass.trainPrecision = myPRClass.trainnpr / float(myPRClass.traindp)
                if myPRClass.trainGood > 0:
                    myPRClass.trainRecall = myPRClass.trainnpr / float(myPRClass.trainGood)
                #print 'trainPrecision: ', myPRClass.trainPrecision
                #print 'trainRecall: ', myPRClass.trainRecall
                if (myPRClass.trainPrecision + myPRClass.trainRecall) > 0:
                    myPRClass.trainF1 = 2 * myPRClass.trainPrecision * myPRClass.trainRecall / (myPRClass.trainPrecision + myPRClass.trainRecall)
                #print 'F1: ', myPRClass.trainF1
                outFile.write('train' + ',' + str(myPRClass.trainPrecision) + ',' + str(myPRClass.trainRecall) \
                              + ',' + str(myPRClass.trainF1) + ',' + str(myPRClass.trainAccuracy)\
                              + ',' + str(myPRClass.T) + ',' + str(myPRClass.k) + ',' \
                              + str(myPRClass.b) + ','\
                              + str(myPRClass.trainMistakes) + ',' + str(t + 1) + '\n')
                logging.info ('train' + ',' + str(myPRClass.trainPrecision) + ',' + str(myPRClass.trainRecall) \
                              + ',' + str(myPRClass.trainF1) + ',' + str(myPRClass.trainAccuracy)\
                              + ',' + str(myPRClass.T) + ',' + str(myPRClass.k) + ',' \
                              + str(myPRClass.b) + ','\
                              + str(myPRClass.trainMistakes) + ',' + str(t + 1) + '\n')

            logging.info('training Complete')

            testTime1 = time.time()
            myPRClass.tp = 0
            myPRClass.fp = 0
            myPRClass.tn = 0
            myPRClass.fn = 0
            myPRClass.testMistakes = 0
            myPRClass.testTotal = 0
            myPRClass.testnpr = 0
            myPRClass.testdp = 0
            myPRClass.testGood = 0
            readGene (thePath, geneTest, kmer, maxReadSize, "testWeight", myPRClass, t)
            testTime = time.time() - testTime1
            outFile.write ('testTime: ' + str(testTime))
            logging.info ('testTime: ' + str(testTime))

            ### report results of testing
            print 'true Positive: ', myPRClass.tp
            print 'false Positive: ', myPRClass.fp
            print 'true Negative: ', myPRClass. tn
            print 'false Negative: ', myPRClass.fn
            if myPRClass.fp + myPRClass.tp != 0:
                print 'positive Success: ', float (myPRClass.tp) / (myPRClass.tp + myPRClass.fp)
            else:
                print 'positive Success: 0'
            if myPRClass.fn + myPRClass.tn != 0:
                print 'negative Success: ', float (myPRClass.tn) / (myPRClass.tn + myPRClass.fn)
            else:
                print 'negative Success: 0'

            myPRClass.testAccuracy = 100 - (100 * (myPRClass.testMistakes / float(myPRClass.testTotal)) )
            #print 'testMistakes = ', myPRClass.testMistakes, '   : test success = ', \
                  #myPRClass.testAccuracy, '%'

            if myPRClass.testdp > 0:
                myPRClass.testPrecision = myPRClass.testnpr / float(myPRClass.testdp)
            myPRClass.testRecall = myPRClass.testnpr / float(myPRClass.testGood)
            #print 'testPrecision: ', myPRClass.testPrecision
            #print 'testRecall: ', myPRClass.testRecall
            if (myPRClass.testPrecision + myPRClass.testRecall) > 0:
                myPRClass.testF1 = 2 * myPRClass.testPrecision * myPRClass.testRecall / (myPRClass.testPrecision + myPRClass.testRecall)

            outFile.write('test' + ',' + str(myPRClass.testPrecision) + ',' + str(myPRClass.testRecall) \
                          + ',' + str(myPRClass.testF1) + ',' + str(myPRClass.testAccuracy)\
                          + ',' + str(myPRClass.T) + ',' + str(myPRClass.k) + ',' \
                          + str(myPRClass.b) + ','\
                          + str(myPRClass.testMistakes) + ',' + str(t + 1) + '\n')
            logging.info('test' + ',' + str(myPRClass.testPrecision) + ',' + str(myPRClass.testRecall) \
                          + ',' + str(myPRClass.testF1) + ',' + str(myPRClass.testAccuracy)\
                          + ',' + str(myPRClass.T) + ',' + str(myPRClass.k) + ',' \
                          + str(myPRClass.b) + ','\
                          + str(myPRClass.testMistakes) + ',' + str(t + 1) + '\n')

    endTime = time.time()
    print 'Perceptron runTime: ', endTime - startTime
    outFile.write (str(endTime - startTime))
    logging.info ('Perceptron runTime: ' + str(endTime - startTime) + '\n')


    ### end test perceptron
    '''
    ### begin pegasos:

    #####
    outFile.write ('pegasos'  + '\n' + '\n')
    outFile.write ('Gene left out: ' + geneTest + '\n')
    logging.info ('Gene left out: ' + geneTest + '\n')
    outFile.write('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
                  'Accuracy' + ',' + 'T' + ',' + 'K' + ',' \
                  + 'Bias' + ','\
                  + 'Mistakes' + ',' + 'Iteration' + '\n')
    logging.info('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
                  'Accuracy' + ',' + 'T' + ',' + 'K' + ',' \
                  + 'Bias' + ','\
                  + 'Mistakes' + ',' + 'Iteration' + '\n')

    print 'running pegasos . . . '
    logging.info ('running pegasos . . . ' + '\n')

    startTime = time.time()

    for t in range(myPGClass.T):

        #print 'iteration: ', t
        #####
        for b in biasVariableList:
            myPGClass.b = b
            for gene in geneList:
                myPGClass.trainTotal = 0
                myPGClass.trainnpr = 0
                myPGClass.traindp = 0
                myPGClass.trainGood = 0
                myPGClass.trainMistakes = 0

                readGene (thePath, gene, kmer, maxReadSize, "pegasos", myPGClass, t)

                ### report results of training

                myPGClass.trainAccuracy = 100 - (100 * (myPGClass.trainMistakes / float(myPGClass.trainTotal)) )
                #print 'trainMistakes = ', myPGClass.trainMistakes, '   : train success = ', \
                      #myPGClass.trainAccuracy, '%'
                if myPGClass.traindp > 0:
                    myPGClass.trainPrecision = myPGClass.trainnpr / float(myPGClass.traindp)
                if myPGClass.trainGood > 0:
                    myPGClass.trainRecall = myPGClass.trainnpr / float(myPGClass.trainGood)
                #print 'trainPrecision: ', myPGClass.trainPrecision
                #print 'trainRecall: ', myPGClass.trainRecall
                if (myPGClass.trainPrecision + myPGClass.trainRecall) > 0:
                    myPGClass.trainF1 = 2 * myPGClass.trainPrecision * myPGClass.trainRecall / (myPGClass.trainPrecision + myPGClass.trainRecall)
                #print 'F1: ', myPGClass.trainF1
                outFile.write('train' + ',' + str(myPGClass.trainPrecision) + ',' + str(myPGClass.trainRecall) \
                              + ',' + str(myPGClass.trainF1) + ',' + str(myPGClass.trainAccuracy)\
                              + ',' + str(myPGClass.T) + ',' + str(myPGClass.k) + ',' \
                              + str(myPGClass.b) + ','\
                              + str(myPGClass.trainMistakes) + ',' + str(t + 1) + '\n')
                logging.info ('train' + ',' + str(myPGClass.trainPrecision) + ',' + str(myPGClass.trainRecall) \
                              + ',' + str(myPGClass.trainF1) + ',' + str(myPGClass.trainAccuracy)\
                              + ',' + str(myPGClass.T) + ',' + str(myPGClass.k) + ',' \
                              + str(myPGClass.b) + ','\
                              + str(myPGClass.trainMistakes) + ',' + str(t + 1) + '\n')

        logging.info('training Complete')

        myPGClass.testTotal = 0
        myPGClass.testnpr = 0
        myPGClass.testdp = 0
        myPGClass.testGood = 0
        myPGClass.testMistakes = 0
        myPGClass.aVal = 0
        myPGClass.bVal = 0
        myPGClass.cVal = 0
        myPGClass.dVal = 0

        readGene (thePath, geneTest, kmer, maxReadSize, "testWeight", myPGClass, t)
        print 'iteration: ', t + 1
        print myPGClass.aVal
        print myPGClass.bVal
        print myPGClass.cVal
        print myPGClass.dVal
        print 'Good rate: ', 100 * (float(myPGClass.aVal) / (myPGClass.aVal + myPGClass.bVal))
        print 'Bad rate: ', 100 * (float(myPGClass.dVal) / (myPGClass.cVal + myPGClass.dVal))
        print 'overall success rate: ', 100 * (float (myPGClass.aVal + myPGClass.dVal) / (myPGClass.aVal + myPGClass.bVal + myPGClass.cVal + myPGClass.dVal))

        ### report results of testing

        myPGClass.testAccuracy = 100 - (100 * (myPGClass.testMistakes / float(myPGClass.testTotal)) )
        #print 'testMistakes = ', myPGClass.testMistakes, '   : test success = ', \
              #myPGClass.testAccuracy, '%'

        if myPGClass.testdp > 0:
            myPGClass.testPrecision = myPGClass.testnpr / float(myPGClass.testdp)
        myPGClass.testRecall = myPGClass.testnpr / float(myPGClass.testGood)
        #print 'testPrecision: ', myPGClass.testPrecision
        #print 'testRecall: ', myPGClass.testRecall
        if (myPGClass.testPrecision + myPGClass.testRecall) > 0:
            myPGClass.testF1 = 2 * myPGClass.testPrecision * myPGClass.testRecall / (myPGClass.testPrecision + myPGClass.testRecall)

        outFile.write('test' + ',' + str(myPGClass.testPrecision) + ',' + str(myPGClass.testRecall) \
                      + ',' + str(myPGClass.testF1) + ',' + str(myPGClass.testAccuracy)\
                      + ',' + str(myPGClass.T) + ',' + str(myPGClass.k) + ',' \
                      + str(myPGClass.b) + ','\
                      + str(myPGClass.testMistakes) + ',' + str(t + 1) + '\n')
        logging.info('test' + ',' + str(myPGClass.testPrecision) + ',' + str(myPGClass.testRecall) \
                      + ',' + str(myPGClass.testF1) + ',' + str(myPGClass.testAccuracy)\
                      + ',' + str(myPGClass.T) + ',' + str(myPGClass.k) + ',' \
                      + str(myPGClass.b) + ','\
                      + str(myPGClass.testMistakes) + ',' + str(t + 1) + '\n')

    endTime = time.time()
    print 'pegasos runTime: ', endTime - startTime
    outFile.write (str(endTime - startTime))
    logging.info ('pegasos runTime: ' + str(endTime - startTime) + '\n')


    ### end test pegasos
    #####
    '''
    '''
    ### Add a header for Pegasos:
    outFile.write ('\n' + '\n' + 'Pegasos' + '\n' + '\n')
    outFile.write('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
                  'Accuracy' + ',' + 'T' + ',' + 'K' + ',' + 'Lamda' + ',' \
                  + 'Mistakes' + ',' + 'Iteration' + '\n')

    print 'running Pegasos . . .'
    logging.info('running Pegasos . . .' + '\n')

    ### record start time to determine runtime of algorithm
    startTime = time.time()

    ### run for each value of T:
    for i in myPGClass.Tlist:
        myPGClass.T = i
        ### run for each value of k:
        for j in myPGClass.klist:
            myPGClass.k = j
            for l in myPGClass.lamList:
                myPGClass.lam = l

                ### repeat for T iterations:
                for t in range(myPGClass.T):
                    #myPGClass.pegasosBatch(trainDataSet, testDataSet, t)

                    outFile.write('train' + ',' + str(myPGClass.trainPrecision) + ',' + str(myPGClass.trainRecall) \
                                  + ',' + str(myPGClass.trainF1) + ',' + str(myPGClass.trainAccuracy)\
                                  + ',' + str(myPGClass.T) + ',' + str(myPGClass.k) + ',' \
                                  + str(myPGClass.lam) + ',' + str(myPGClass.trainMistakes) + ',' + str(t + 1) + '\n')

                    #myPGClass.testWeight(testDataSet)
                    outFile.write('test' + ',' + str(myPGClass.testPrecision) + ',' + str(myPGClass.testRecall) \
                                  + ',' + str(myPGClass.testF1) + ',' + str(myPGClass.testAccuracy)\
                                  + ',' + str(myPGClass.T) + ',' + str(myPGClass.k) + ',' \
                                  + str(myPGClass.lam) + ',' + str(myPGClass.testMistakes) + '\n')


    endTime = time.time()
    print 'Pegasos runTime: ', endTime - startTime
    outFile.write (str(endTime - startTime))
    logging.info('Pegasos runTime: ' + str(endTime - startTime) + '\n')

    ### end pegasos
    '''

    print 'closing output file . . .'
    logging.info ('closing output file . . .' + '\n')

    outFile.close()

### end individualProcess function

### Run each N in parallel:

processes = [mp.Process(target=individualProcess, args=(k, l)) for k in kmerList for l in subList]

for p in processes:
    p.start()

for p in processes:
    p.join()

print '*** program complete ***'
logging.info ('*** program complete ***' + '\n')
