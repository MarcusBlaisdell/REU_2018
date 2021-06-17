###
# count unique indexes
###

### loadDatafunction

word = ''
value = ''
line = []
dataList = []
uniqueDict = {}
label = 0
readFile = open("/home/marcus/Data/Genomic-Data/testing_Set-3", "r")

stringOne = readFile.readline()

if stringOne[0] == 'f':
    ### parse out header from first line by looking for occurance of specific characters:
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

### now, parse data into a list:
while string:
    i = 0

    if string[len(string) - 3] == '-':
        label = -1
    else:
        label = 1

    ### value comes after the colon,
    ### index is before
    while i < len(string):
        if string[i] == ':':
            i += 1
            while string[i] != ' ':
                value += string[i]
                i += 1
            tuple = [int(word), int(value)]

            ### add unique indexes to the dictionary

            if uniqueDict.get(int(word),'--') == '--':
                uniqueDict[int(word)] = 1

            line.append(tuple)
            word = ''
            value = ''

        else:
            word += string[i]
        i += 1

    lineTuple = [line, label]
    #self.trainData.append(lineTuple)
    dataList.append(lineTuple)
    line = []
    word = ''
    string = readFile.readline()

readFile.close()

### end loadData functions

writeFile = open("uniqueIndexes", "w")

keyList = uniqueDict.keys()
#keyList.sort()
print 'number of unique indexes: ', len(keyList)
print 'largest index: ', max(keyList)

for key in keyList:
    writeFile.write(str(key) + '\n')

writeFile.close()
