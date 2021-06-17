##########################################
### PEGASOS Algorithm Test
### for genomic data
### Marcus Blaisdell
##########################################

from p_Class import *

myClass = p_Class()

myClass.loadData("/home/marcus/Data/Genomic-Data/train-11", myClass.trainData)
#myClass.loadData("/home/marcus/Data/Genomic-Data/shortEDL-validation.txt", myClass.validationData)
myClass.loadData("/home/marcus/Data/Genomic-Data/test-11", myClass.testData)
myClass.countGood()

### open the results file:
#outFile = open("/home/marcus/Documents/LinuxShare/Results/Pegasos/results_variables-T.csv", "w")
outFile = open("results_variables-T.csv", "w")

### Add a header:
outFile.write('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
              'Accuracy' + ',' + 'T' + ',' + 'K' + ',' + 'Lamda' + ',' \
              + 'Mistakes' + ',' + 'Iteration' + '\n')

### record start time to determine runtime of algorithm
startTime = time.time()

### run for each value of T:
for i in myClass.Tlist:
    myClass.T = i
    ### run for each value of k:
    for j in myClass.klist:
        myClass.k = j
        for l in myClass.lamList:
            myClass.lam = l

            '''
            # initialize weight vector, w to zero:

            kmerSize = ((4**myClass.kmer) + 1)
            myClass.w = np.zeros(kmerSize, dtype=float)
            #myClass.w = np.zeros(10000000, dtype = float)
            '''

            ### repeat for T iterations:
            for t in range(myClass.T):
                myClass.pegasosBatchBiasTwo(t)

                outFile.write('train' + ',' + str(myClass.trainPrecision) + ',' + str(myClass.trainRecall) \
                              + ',' + str(myClass.trainF1) + ',' + str(myClass.trainAccuracy)\
                              + ',' + str(myClass.T) + ',' + str(myClass.k) + ',' \
                              + str(myClass.lam) + ',' + str(myClass.trainMistakes) + ',' + str(t + 1) + '\n')

                myClass.testWeight()
                outFile.write('test' + ',' + str(myClass.testPrecision) + ',' + str(myClass.testRecall) \
                              + ',' + str(myClass.testF1) + ',' + str(myClass.testAccuracy)\
                              + ',' + str(myClass.T) + ',' + str(myClass.k) + ',' \
                              + str(myClass.lam) + ',' + str(myClass.testMistakes) + '\n')

'''
for i in myClass.Tlist:
    myClass.T = i
    for j in myClass.lamList:
        myClass.lam = j
        for k in myClass.klist:
            myClass.k = k

            # repeat test for x cycles:
            x = 1
            for m in range (x):
                #pegasosStart = time.time()
                myClass.nonZeroDict.clear()
                myClass.pegasosNew()

                outFile.write('train' + ',' + str(myClass.trainPrecision) + ',' + str(myClass.trainRecall) \
                              + ',' + str(myClass.trainF1) + ',' + str(myClass.trainAccuracy)\
                              + ',' + str(myClass.T) + ',' + str(myClass.k) + ',' \
                              + str(myClass.lam) + ',' + str(myClass.trainMistakes) + '\n')

                #pegasosEnd = time.time()
                #print 'Pegasos: ', pegasosEnd - pegasosStart

                #validationStart = time.time()

                myClass.validationWeight()
                outFile.write('validation' + ',' + str(myClass.validationPrecision) + ',' + str(myClass.validationRecall) \
                              + ',' + str(myClass.validationF1) + ',' + str(myClass.validationAccuracy)\
                              + ',' + str(myClass.T) + ',' + str(myClass.k) + ',' \
                              + str(myClass.lam) + ',' + str(myClass.validationMistakes) + '\n')

                #validationEnd = time.time()
                #print 'validation: ', validationEnd - validationStart

                #testStart = time.time()

                myClass.testWeight()
                outFile.write('test' + ',' + str(myClass.testPrecision) + ',' + str(myClass.testRecall) \
                              + ',' + str(myClass.testF1) + ',' + str(myClass.testAccuracy)\
                              + ',' + str(myClass.T) + ',' + str(myClass.k) + ',' \
                              + str(myClass.lam) + ',' + str(myClass.testMistakes) + '\n')

                #testEnd = time.time()
                #print 'test: ', testEnd - testStart
'''
endTime = time.time()
print 'runTime: ', endTime - startTime
outFile.write (str(endTime - startTime) + '\n')
#print 'bias: ', myClass.b

outFile.close()

#timeArray.append(endTime - startTime)
'''
totalTime = 0.0
for x in timeArray:
    totalTime += x
avgTime = totalTime / 10
print 'average: ' , totalTime / 10
#avgFile = open("time.txt", "w")
#avgFile.write ('average: ' + str(avgTime) + '\n')
#avgFile.close()
'''
