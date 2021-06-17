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

### sign function
def sign(number):
    if number > 0:
        return 1
    else:
        return -1

### end sign function
