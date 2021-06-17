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
                        i += 1
                except:
                    i += 1

            string = stringOne[i:]
            ### end parse Header
        else:
            string = stringOne

        ### now, parse data into a list

        ### use i to iterate through the string:

        i = 0

        ### The label is at the end:

        if (len(string)) > 2:
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
        else:
            return 0

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

def readGene (thePath, gene, kmer, maxReadSize, functionName, myClass, t, logging, writeFile):
    logging.info ('readGene called' + '\n')
    fileName = thePath + gene + '-' + str(kmer)
    #print fileName
    readFile = open (fileName, 'r')
    readFlag = 1
    while readFlag == 1:
        runningList = []
        readFlag = readData (runningList, maxReadSize, readFile)
        if functionName == 'perceptron':
            myClass.perceptron (runningList)
        if functionName == 'testWeight':
            myClass.testWeight (runningList, t)
        if functionName == 'testWeightPair':
            myClass.testWeightPair (runningList, t)
        if functionName == 'newTest':
            myClass.newTest (runningList, t, writeFile)

    readFile.close()

### End of file
