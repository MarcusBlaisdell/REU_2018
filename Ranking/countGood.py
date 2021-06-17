geneList = ['0157', 'Co6114', 'CH611', '2016C']

thePath = "/Users/MarcusBlaisdell/Documents/LinuxShare/tenK/"
for gene in geneList:
    good = 0
    bad = 0
    theFile = thePath + gene + "-11"
    readFile = open(theFile, "r")

    string = readFile.readline()

    while string:
        if len(string) > 5:
            label = int(string.split(" ")[-1])
            if label == 1:
                good += 1
            else:
                bad += 1

        string = readFile.readline()

    print gene, ': good: ', good, ' : bad: ', bad
    readFile.close()
