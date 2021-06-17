##########################################
### Common functions
###
### Marcus Blaisdell
##########################################

### loadDatafunction

def loadData(fileName, dataName):
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
                tuple = [int(word), int(value)]
                line.append(tuple)
                word = ''
                value = ''

            else:
                word += string[i]
            i += 1

        lineTuple = [line, label]

        dataName.append(lineTuple)
        line = []
        word = ''
        string = readFile.readline()

    readFile.close()

### end loadData functions

#####
### loadDataOne function

def loadDataOne (fileName):
    word = ''
    value = ''
    line = []
    label = 0
    returnArray = []

    stringOne = fileName.readline()

    if stringOne:
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
        ### parse data into an array:
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
                tuple = [int(word), int(value)]
                line.append(tuple)
                word = ''
                value = ''

            else:
                word += string[i]
            i += 1

        returnArray = [line, label]

        line = []
        word = ''

        return returnArray

### end loadDataOne function
#####

### sign function
def sign(number):
    if number > 0:
        return 1
    else:
        return -1

### end sign function

### Process a single line read from a file:

def processLine(readFile, dataName):
    word = ''
    value = ''
    line = []
    label = 0

    ### read a line from the file:
    stringOne = readFile.readline()

    ### if the string is empty, return 0,
    ### otherwise, process the string

    if not stringOne:
        return 0
    else:
        ### if this is the first line of the file,
        ### we need to parse out the header:

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

        ### use i to iterate through the string:

        i = 0

        ### The label is at the end:

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
                tuple = [int(word), int(value)]
                line.append(tuple)
                word = ''
                value = ''

            else:
                word += string[i]
            i += 1

        lineTuple = [line, label]

        #dataName.append(lineTuple)
        line = []
        word = ''
        #string = readFile.readline()

    return lineTuple

### Read data from file in maxReadSize batches:

def readData (dataName, maxReadSize, readFile):
    for x in range(maxReadSize):
        #string = readFile.readline()
        string = processLine (readFile, dataName)

        if not string:
            return 0
        else:
            dataName.append (string)

    ### if read was good, return good
    return 1

### Read one gene and process it through the requested function:

def readGene (thePath, gene, kmer, maxReadSize, functionName, myClass, t):
    fileName = thePath + gene + '-' + str(kmer)
    # fileName
    '''
    #####
    readFile = open ('/Users/MarcusBlaisdell/Documents/LinuxShare/tenK/0157-11', 'r')
    testList = []
    rFlag = readData (testList, maxReadSize, readFile)
    readFile.close()
    #####
    '''
    readFile = open (fileName, 'r')
    readFlag = 1
    while readFlag == 1:
        runningList = []
        readFlag = readData (runningList, maxReadSize, readFile)
        if functionName == 'perceptron':
            myClass.perceptron (runningList)
        if functionName == 'pegasos':
            myClass.pegasosBatch (runningList, t)
            '''
            myClass.testWeight (testList, t)
            ### report results of testing

            outFile.write ('pegasos'  + '\n' + '\n')
    outFile.write ('Gene left out: ' + geneTest + '\n')
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
    for biasVariable in biasVariableList:
        myClass.b = biasVariable
        #####
        for t in range(myPGClass.T):
            #print 'iteration: ', t
            #####
            for gene in geneList:
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

            readGene (thePath, geneTest, kmer, maxReadSize, "testWeight", myPGClass, t)
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


    ### end test pegasosoutFile.write ('pegasos'  + '\n' + '\n')
    outFile.write ('Gene left out: ' + geneTest + '\n')
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
    for biasVariable in biasVariableList:
        myPGClass.b = biasVariable
        #####
        for t in range(myPGClass.T):
            #print 'iteration: ', t
            #####
            for gene in geneList:
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

            readGene (thePath, geneTest, kmer, maxReadSize, "testWeight", myPGClass, t)
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


    ### end test pegasos.testAccuracy = 100 - (100 * (myPRClass.testMistakes / float(myPRClass.testTotal)) )
            #print 'testMistakes = ', myPRClass.testMistakes, '   : test success = ', \
                  #myPRClass.testAccuracy, '%'

            if myPRClass.testdp > 0:
                myPRClass.testPrecision = myPRClass.testnpr / float(myPRClass.testdp)
            myPRClass.testRecall = myPRClass.testnpr / float(myPRClass.testGood)
            #print 'testPrecision: ', myPRClass.testPrecision
            #print 'testRecall: ', myPRClass.testRecall
            if (myPRClass.testPrecision + myPRClass.testRecall) > 0:
                myPRClass.testF1 = 2 * myPRClass.testPrecision * myPRClass.testRecall / (myPRClass.testPrecision + myPRClass.testRecall)

            print ('test' , ',' , str(myPRClass.testPrecision) , ',' ,  str(myPRClass.testRecall) \
                          , ',' , str(myPRClass.testF1) , ',' , str(myPRClass.testAccuracy)\
                          , ',' , str(myPRClass.T) , ',' , str(myPRClass.k) , ',' \
                          , str(myPRClass.b) , ','\
                          , str(myPRClass.testMistakes) , ',' , str(t + 1) , '\n')
        '''
        if functionName == 'testWeight':
            myClass.testWeight (runningList, t)

    readFile.close()

### End of file
