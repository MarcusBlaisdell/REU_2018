##########################################
### Perceptron Test
### Marcus Blaisdell
##########################################

from p_Class import *
#import time

myClass = p_Class()

myClass.loadData("/Users/MarcusBlaisdell/Documents/LinuxShare/tenK/2016C-11", myClass.trainData)
myClass.loadData("/Users/MarcusBlaisdell/Documents/LinuxShare/tenK/0157-11", myClass.testData)
'''
myClass.loadData("/home/marcus/Data/Genomic-Data/training_Set-1", myClass.trainData)
#myClass.loadData("/home/marcus/Data/Genomic-Data/shortEDL-test.txt", myClass.validationData)
myClass.loadData("/home/marcus/Data/Genomic-Data/testing_Set-1", myClass.testData)
'''
'''
myClass.loadData("/home/marcus/Data/Genomic-Data/training", myClass.trainData)
myClass.loadData("/home/marcus/Data/Genomic-Data/testing", myClass.testData)
'''
myClass.countGood()

#outFile = open("/home/marcus/Documents/LinuxShare/Results/Perceptron/results_Set-1-no-bias.csv", "w")
outFile = open("results_Set-1-no-bias.csv", "w")

outFile.write('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
              'Accuracy' + ',' + 'T' + ',' + 'K' + ',' \
              + 'Mistakes' + ',' + 'Iteration' + '\n')

startTime = time.time()
myClass.perceptronBias(outFile)
endTime = time.time()
print 'runTime: ', endTime - startTime

outFile.close()
