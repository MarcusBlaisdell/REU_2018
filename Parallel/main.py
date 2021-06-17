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

### create array of the values to label the training, testing, results files:

#kmerList = [11, 13, 15, 17]
kmerList = [11]
#subList = ['a', 'b', 'c', 'd']
subList = ['a']
#biasVariableList = [1, 10, 50, 100]
biasVariableList = [0]

### run the comparison for each training/testing set:

thePath = "/home/marcus/Data/Genomic-Data/Testing/0157-1x/"
#thePath = "/data/doppa/users/mblaisdell/genome-assembly-datasets/Testing/"
#resultsStringPre = "/home/marcus/Documents/LinuxShare/Results/resultszU-A-b-1-100_"
resultsStringPre = "results_U-A_L2-0157-1x_"
resultsStringPost = ".csv"

### create empty status file:
statusFile = open("status.txt", "w")
statusFile.close()

def individualProcess(kmer, letter):
    myPGClass = pg_Class()
    myPRClass = pr_Class()

    ### The datasets are the same for each function, create them once to save memmory:
    trainDataSet = []
    testDataSet = []

    myPGClass.kmer = kmer
    myPRClass.kmer = kmer

    trainData = thePath + str(kmer) + letter + '-train'
    testData = thePath + str(kmer) + letter + '-test'
    resultsFile = resultsStringPre + str(kmer) + letter + resultsStringPost
    print 'resultsFile: ', resultsFile

    statusFile = open ("status.txt", "a")
    print 'set: ', str(kmer) + letter
    statusFile.write('set: ' + str(str(kmer) + letter) + '\n')

    print 'loading data:'
    statusFile.write ('loading data:' + '\n')

    trainDataSet = []
    testDataSet = []

    loadData(trainData, trainDataSet)
    loadData(testData, testDataSet)

    print 'load complete'
    statusFile.write ('load complete' + '\n')

    print 'opening results file:'
    statusFile.write ('opening results file:' + '\n')
    outFile = open(resultsFile, "w")

    myPRClass.countGood(trainDataSet, testDataSet)

    outFile.write ('Perceptron'  + '\n' + '\n')
    outFile.write('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
                  'Accuracy' + ',' + 'T' + ',' + 'K' + ',' \
                  + 'Bias' + ','\
                  + 'Mistakes' + ',' + 'Iteration' + '\n')

    print 'running perceptron . . . '
    statusFile.write ('running perceptron . . . ' + '\n')
    statusFile.close()

    startTime = time.time()
    for biasVariable in biasVariableList:
        myPRClass.b = biasVariable
        myPRClass.perceptron(trainDataSet, testDataSet, outFile)
    endTime = time.time()
    print 'Perceptron runTime: ', endTime - startTime
    outFile.write (str(endTime - startTime))
    statusFile = open("status.txt", "a")
    statusFile.write('Perceptron runTime: ' + str(endTime - startTime) + '\n')

    ### end test perceptron

    '''
    myPGClass.countGood(trainDataSet, testDataSet)

    ### Add a header for Pegasos:
    outFile.write ('\n' + '\n' + 'Pegasos' + '\n' + '\n')
    outFile.write('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
                  'Accuracy' + ',' + 'T' + ',' + 'K' + ',' + 'Lamda' + ',' \
                  + 'Mistakes' + ',' + 'Iteration' + '\n')

    print 'running Pegasos . . .'
    statusFile.write('running Pegasos . . .' + '\n')
    statusFile.close()

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
    statusFile = open("status.txt", "a")
    statusFile.write('Pegasos runTime: ' + str(endTime - startTime) + '\n')
    '''

    print 'closing output file . . .'
    statusFile.write('closing output file . . .' + '\n')
    statusFile.close()

    outFile.close()

    #print 'zeroCount: ', myPRClass.zeroCount ### Testing

### end individualProcess function

processes = [mp.Process(target=individualProcess, args=(k, l)) for k in kmerList for l in subList]

for p in processes:
    p.start()

for p in processes:
    p.join()


'''
for kmer in kmerList:
    myPGClass = pg_Class()
    myPRClass = pr_Class()

    ### The datasets are the same for each function, create them once to save memmory:
    trainDataSet = []
    testDataSet = []

    myPGClass.kmer = kmer
    myPRClass.kmer = kmer
    for letter in subList:
        trainData = thePath + str(kmer) + letter + '-train'
        testData = thePath + str(kmer) + letter + '-test'
        resultsFile = resultsStringPre + str(kmer) + letter + resultsStringPost
        print 'resultsFile: ', resultsFile

        statusFile = open ("status.txt", "a")
        print 'set: ', str(kmer) + letter
        statusFile.write('set: ' + str(str(kmer) + letter) + '\n')

        print 'loading data:'
        statusFile.write ('loading data:' + '\n')

        trainDataSet = []
        testDataSet = []

        loadData(trainData, trainDataSet)
        loadData(testData, testDataSet)
        print 'load complete'
        statusFile.write ('load complete' + '\n')

        print 'opening results file:'
        statusFile.write ('opening results file:' + '\n')
        outFile = open(resultsFile, "w")

        myPRClass.countGood(trainDataSet, testDataSet)

        outFile.write ('Perceptron'  + '\n' + '\n')
        outFile.write('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
                      'Accuracy' + ',' + 'T' + ',' + 'K' + ',' \
                      + 'Mistakes' + ',' + 'Iteration' + '\n')

        print 'running perceptron . . . '
        statusFile.write ('running perceptron . . . ' + '\n')
        statusFile.close()

        startTime = time.time()
        myPRClass.perceptron(trainDataSet, testDataSet, outFile)
        endTime = time.time()
        print 'Perceptron runTime: ', endTime - startTime
        outFile.write (str(endTime - startTime))
        statusFile = open("status.txt", "a")
        statusFile.write('Perceptron runTime: ' + str(endTime - startTime) + '\n')

        ### end test perceptron

        myPGClass.countGood(trainDataSet, testDataSet)


        ### Add a header for Pegasos:
        outFile.write ('\n' + '\n' + 'Pegasos' + '\n' + '\n')
        outFile.write('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
                      'Accuracy' + ',' + 'T' + ',' + 'K' + ',' + 'Lamda' + ',' \
                      + 'Mistakes' + ',' + 'Iteration' + '\n')

        print 'running Pegasos . . .'
        statusFile.write('running Pegasos . . .' + '\n')
        statusFile.close()

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
                        myPGClass.pegasosBatch(trainDataSet, testDataSet, t)

                        outFile.write('train' + ',' + str(myPGClass.trainPrecision) + ',' + str(myPGClass.trainRecall) \
                                      + ',' + str(myPGClass.trainF1) + ',' + str(myPGClass.trainAccuracy)\
                                      + ',' + str(myPGClass.T) + ',' + str(myPGClass.k) + ',' \
                                      + str(myPGClass.lam) + ',' + str(myPGClass.trainMistakes) + ',' + str(t + 1) + '\n')

                        myPGClass.testWeight(trainDataSet, testDataSet)
                        outFile.write('test' + ',' + str(myPGClass.testPrecision) + ',' + str(myPGClass.testRecall) \
                                      + ',' + str(myPGClass.testF1) + ',' + str(myPGClass.testAccuracy)\
                                      + ',' + str(myPGClass.T) + ',' + str(myPGClass.k) + ',' \
                                      + str(myPGClass.lam) + ',' + str(myPGClass.testMistakes) + '\n')


        endTime = time.time()
        print 'Pegasos runTime: ', endTime - startTime
        outFile.write (str(endTime - startTime))
        statusFile = open("status.txt", "a")
        statusFile.write('Pegasos runTime: ' + str(endTime - startTime) + '\n')

        print 'closing output file . . .'
        statusFile.write('closing output file . . .' + '\n')
        statusFile.close()

        outFile.close()
'''

print '*** program complete ***'
statusFile = open("status.txt", "a")
statusFile.write('*** program complete ***' + '\n')
statusFile.close()
