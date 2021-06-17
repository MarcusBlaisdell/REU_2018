##########################################
### ranking Test
### for genomic data
### Marcus Blaisdell
### Parallel Processing
###
##########################################

from sklearn.metrics import precision_score
import multiprocessing as mp
from pg_Class import *
from pr_Class import *
from r_Class import *
from functions import *
import numpy as np
import logging
import sys
import time

### create array of the values to label the training, testing, results files:

#geneList = ['0157', '2016C', 'CH611', 'Co6114', 'ED1a', 'EDL933-1', 'FAP1', '_isolate102', 'RS76', 'UMN026']
#geneList = ['2016C', 'CH611', 'Co6114']
geneList = ['0157', 'CH611', '2016C']
#geneList = ['2016C']
geneTest = 'Co6114'
#geneTest = '2016C'
#kmerList = [11, 13, 15, 17]
kmerList = [11]
#subList = ['a', 'b', 'c', 'd']
subList = ['a']
#biasVariableList = [1, 10, 50, 100]
biasVariableList = [0]

maxReadSize = 500000
sampleSizeK = 10000

### run the comparison for each training/testing set:

#thePath = "/home/marcus/Data/Genomic-Data/Testing/"
#thePath = "/data/doppa/users/mblaisdell/genome-assembly-datasets/"
thePath = "/Users/MarcusBlaisdell/Documents/LinuxShare/tenK/"
#thePath = ""
#resultsStringPre = "/home/marcus/Documents/LinuxShare/Results/resultszU-A-b-1-100_"
resultsStringPre = "results_ranking_11_Tau-20-Point8-b_" + geneTest + "_"
resultsStringPost = ".csv"

### create empty log file:
logging.basicConfig(filename='run.log-11_Tau-20-Point8-b',level=logging.DEBUG)
#logging.info("Log File created")
logging.info (geneTest + '\n')

def individualProcess(kmer, letter):
    myPGClass = pg_Class ()
    myPRClass = pr_Class ()
    myRClass = r_Class ()

    myPGClass.kmer = kmer
    myRClass.kmer = kmer
    myRClass.k = sampleSizeK

    resultsFile = resultsStringPre + str(kmer) + resultsStringPost
    print 'resultsFile: ', resultsFile
    #logging.info ('results File: ' + resultsFile)

    print 'set: ', str(kmer)
    logging.info ('set: ' + str(str(kmer)) + '\n')

    print 'loading data:'
    #logging.info('loading data')

    print 'opening results file:'
    #logging.info ('opening results file:' + '\n')
    outFile = open(resultsFile, "w")
    writeFile = open('evaluation_Tau-20-Point8.csv', 'w')
    writeFile.write ('tau' + ',' + 'threshold' + ',' + 'target precision' + ',' + 'tpfb' + ','\
                    + 'fpfb' + ',' + 'yPlus' + ',' + 'precision' + ',' + 'recall' + '\n')

    outFile.write ('ranking'  + '\n' + '\n')
    outFile.write ('Gene left out: ' + geneTest + '\n')
    outFile.write('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
                  'Accuracy' + ',' + 'T' + ',' + 'K' + ',' \
                  + 'Bias' + ','\
                  + 'Mistakes' + ',' + 'Iteration' + '\n')
    logging.info('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
                  'Accuracy' + ',' + 'T' + ',' + 'K' + ',' \
                  + 'Bias' + ','\
                  + 'Mistakes' + ',' + 'Iteration' + '\n')

    print 'running ranking . . . '
    #logging.info ('running ranking . . . ' + '\n')

    startTime = time.time()
    for biasVariable in biasVariableList:
        myRClass.b = biasVariable
        #####

        for t in range(myRClass.T):
            writeFile.write('iteration: ' + str(t) + '\n')
            #print 'iteration: ', t
            #####
            #myRClass.scores = []
            #myRClass.labels = []
            for gene in geneList:
                myRClass.fileLength = myRClass.fileList.get(gene)
                trainStart = time.time()
                #myRClass.pairwise (sampleSizeK, thePath, gene, kmer)
                print 'gene: ', gene
                #myRClass.newPairwise (sampleSizeK, thePath, gene, kmer, writeFile, logging)
                myRClass.newPairwise (sampleSizeK, thePath, gene, kmer, writeFile, logging)
                trainEnd = time.time()

                ### report results of training

                myRClass.trainAccuracy = 100 - (100 * (myRClass.trainMistakes / float(myRClass.trainTotal)) )
                if myRClass.traindp > 0:
                    myRClass.trainPrecision = myRClass.trainnpr / float(myRClass.traindp)
                if myRClass.trainGood > 0:
                    #myRClass.trainRecall = myRClass.trainnpr / float(myRClass.trainGood)
                    myRClass.trainRecall = myRClass.trainnpr / float(myRClass.traingoodStar)
                #print 'trainPrecision: ', myRClass.trainPrecision
                #print 'trainRecall: ', myRClass.trainRecall
                if (myRClass.trainPrecision + myRClass.trainRecall) > 0:
                    myRClass.trainF1 = 2 * myRClass.trainPrecision * myRClass.trainRecall / (myRClass.trainPrecision + myRClass.trainRecall)
                #print 'F1: ', myRClass.trainF1
                outFile.write('train' + ',' + str(myRClass.trainPrecision) + ',' + str(myRClass.trainRecall) \
                              + ',' + str(myRClass.trainF1) + ',' + str(myRClass.trainAccuracy)\
                              + ',' + str(myRClass.T) + ',' + str(myRClass.k) + ',' \
                              + str(myRClass.b) + ','\
                              + str(myRClass.trainMistakes) + ',' + str(t + 1) + '\n')
                logging.info ('train' + ',' + str(myRClass.trainPrecision) + ',' + str(myRClass.trainRecall) \
                              + ',' + str(myRClass.trainF1) + ',' + str(myRClass.trainAccuracy)\
                              + ',' + str(myRClass.T) + ',' + str(myRClass.k) + ',' \
                              + str(myRClass.b) + ','\
                              + str(myRClass.trainMistakes) + ',' + str(t + 1) + '\n')
                outFile.write ('train Time: ' + str(trainEnd - trainStart) + '\n')
                logging.info ('train Time: ' + str(trainEnd - trainStart) + '\n')

            #logging.info('training Complete')

            myRClass.saveWeight()


            endTime = time.time()
            print 'Ranking runTime: ', endTime - startTime
            outFile.write (str(endTime - startTime))
            logging.info ('Ranking runTime: ' + str(endTime - startTime) + '\n')

            ### test ranking

            fileSize = myRClass.fileList.get(geneTest)
            testTime1 = time.time()
            myRClass.newTest (sampleSizeK, thePath, geneTest, kmer, writeFile, logging)
            testTime = time.time() - testTime1
            outFile.write (str(testTime) + '\n')
            logging.info ('test time: ' + str(testTime) + '\n')

            '''
            scores = np.array(myRClass.scores)
            labels = np.array(myRClass.labels)
            positive_scores = scores[labels == 1]
            print 'len(positive_scores): ', len(positive_scores)
            threshold = np.percentile (positive_scores, myRClass.target * 100)
            predicted = scores >= threshold
            print 'len(predicted): ', len(predicted)
            print 'predicted: ', predicted
            print 'labels: ', labels
            pScore = precision_score (labels, predicted)
            print 'precision score at ', myRClass.target, ': ', pScore
            outFile.write ('precision score at ' + str(myRClass.target) + ' recall: ' + str(pScore) + '\n')
            logging.info ('precision score at ' + str(myRClass.target) + ' recall: ' + str(pScore) + '\n')
            '''

            ### report results of testing

            myRClass.testAccuracy = 100 - (100 * (myRClass.testMistakes / float(myRClass.testTotal)) )
            #print 'testMistakes = ', myRClass.testMistakes, '   : test success = ', \
                  #myRClass.testAccuracy, '%'
            if myRClass.testdp > 0:
                myRClass.testPrecision = myRClass.testnpr / float(myRClass.testdp)
            if myRClass.testGood > 0:
                myRClass.testRecall = myRClass.testnpr / float(myRClass.testGood)
            #print 'testPrecision: ', myRClass.testPrecision
            #print 'testRecall: ', myRClass.testRecall
            if (myRClass.testPrecision + myRClass.testRecall) > 0:
                myRClass.testF1 = 2 * myRClass.testPrecision * myRClass.testRecall / (myRClass.testPrecision + myRClass.testRecall)
            #print 'F1: ', myRClass.testF1
            outFile.write('test' + ',' + str(myRClass.testPrecision) + ',' + str(myRClass.testRecall) \
                          + ',' + str(myRClass.testF1) + ',' + str(myRClass.testAccuracy)\
                          + ',' + str(myRClass.T) + ',' + str(myRClass.k) + ',' \
                          + str(myRClass.testMistakes) + ',' + str(t + 1) + '\n')
            logging.info ('test' + ',' + str(myRClass.testPrecision) + ',' + str(myRClass.testRecall) \
                          + ',' + str(myRClass.testF1) + ',' + str(myRClass.testAccuracy)\
                          + ',' + str(myRClass.T) + ',' + str(myRClass.k) + ',' \
                          + str(myRClass.testMistakes) + ',' + str(t + 1) + '\n')
            ### end test ranking

            print 'closing output file . . .'
            #logging.info ('closing output file . . .' + '\n')

    outFile.close()
    writeFile.close()

### end individualProcess function

### Run each N in parallel:

processes = [mp.Process(target=individualProcess, args=(k, l)) for k in kmerList for l in subList]

for p in processes:
    p.start()

for p in processes:
    p.join()

print '*** program complete ***'
logging.info ('*** program complete ***' + '\n')
