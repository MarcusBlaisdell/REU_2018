##########################################
### ranking Test
### class
### Marcus Blaisdell
##########################################

import numpy as np
import random
from functions import *

### Create class pw_Class for ranking:
class r_Class ():
    w = {}
    scoreList = []
    eta = 1.0 # eta is learning rate
    T = 3
    k = 1
    b = 0.0 # default biasVariable
    trainMistakes = 0
    trainTotal = 0
    trainGood = 0
    trainnpr = 0 # numerator for trainPrecision and trainRecall
    traindp = 0 # denominator for trainPrecision
    trainPrecision = 0.0
    trainRecall = 0.0
    trainF1 = 0.0
    testMistakes = 0
    testTotal = 0
    testGood = 0
    testnpr = 0 # numerator for trainPrecision and trainRecall
    testdp = 0 # denominator for trainPrecision
    testPrecision = 0.0
    testRecall = 0.0
    testF1 = 0.0
    TruePositive = 0
    FalsePositive = 0
    TrueNegative = 0
    FalseNegative = 0
    underSample = 0
    scores = []
    labels = []
    target = 0.8
    targetList = [0.98, 0.8, 0.5, 0.25, 0.1]
    threshold = 0.0
    tau = 20 # the percent of the top tau
    #tauList = [10, 20, 50, 75, 100]
    tauList = [99]
    traingoodStar = 0
    testgoodStar = 0


    #fileLength = 38000000 # temporary placeholder, replace by actual file size (computed)
    fileLength = 10000
    fileList = {"0157":500000, "2016C":500000, "CH611":500000, "Co6114":500000}
    #fileList = {"0157":39530625, "2016C":37059519, "CH611":34837752, "Co6114":35684248, "ED1a":38008400, "EDL933-1":38957663, "FAP1":35662001, "_isolate102":33888563, "RS76":34293140, "UMN026":38092060}


    def __init__(self):
        pass

    ### pairwise function:

    def pairwise (self, sampleSizeK, thePath, gene, kmer):
        fileName = thePath + gene + '-' + str(kmer)

        # 1. sample K, k-mers:
        # 2. score k-mers using current weight
        # 3. sort them according to scores
        # 4. identify mistaken pairs (good, bad)
        # update weight: w_t+1 = w_t + eta * (x_good - x_bad)

        ### 1. Sample k k-mers:

        ### create a list of size sampleSizeK of random line numbers
        ### from the file:
        sampleList = self.getRandomList(sampleSizeK)
        self.trainTotal = sampleSizeK
        self.trainMistakes = 0

        ### 2. Score k-mers using current weights

        ### read the selected samples from the file into a list:
        ### and calculate their score:

        readFile = open (fileName, 'r')

        r = 0
        self.scoreList = []

        for s in sampleList:
            t = s - r - 1
            ### iterate to the line before the one in s:
            for x in range(t):
                tempString = readFile.readline()
            ### read the sample:
            sampleX = processLine(readFile, 'bogus')
            if sampleX == 0:
                continue
            theScore = self.calcScore (sampleX)
            theLabel = sampleX[1]
            self.scores.append(theScore)
            if theLabel == 1:
                self.labels.append(theLabel)
            else:
                self.labels.append(0)

            scoreTuple = (theScore, sampleX)
            #TODO **replace with insertSort
            self.scoreList.append(scoreTuple)
            r = s

        readFile.close()

        ### 3. sort them according to scores

        ### sort the scores: (This could be replaced by insertSort)
        ### Put the highest scores at the top
        ### to rank them from highest to lowest
        self.scoreList.sort(reverse=True)

        ### 4. identify mistaken pairs (good, bad)

        ### For the top 20% of scores,
        ### check each pair, if we have an unmatched pair, (good,bad), (bad,good),
        ### update the weight with their difference:

        ### get the total number of good from the list of samples: (for metrics)

        for k in range(len(self.scoreList)):
            if self.scoreList[k][1][1] == 1:
                self.trainGood += 1
                self.traingoodStar += 1

        ### a mistake is an irrelevant record scored higher than a relevant record:
        ### traverse the list and count irrelevant records until a relevant record
        ### is located then add that to the mistakes sum:
        ###
        ### *** New ***, as of August 17, the number of mistakes is the count of
        ### irrelevant records in the top tau:
        for k in range(int (len(self.scoreList) * self.tau / 100.0)):
            if self.scoreList[k][1][1] == 1:
                self.trainMistakes += 1
        '''
        k = 0
        kFlag = 0

        while kFlag == 0:
            count = 1
            if self.scoreList[k][1][1] == -1:
                l = k + 1
                lFlag = 0
                while lFlag == 0:
                #while (l < len(self.scoreList)):
                    if (self.scoreList[l][1][1] == 1):
                        #self.updateWeight (self.scoreList[k], self.scoreList[l])
                        self.trainMistakes += count
                        k = l + 1 # only highest irrelevant
                        ### if we reach the end of the list, exit
                        if k >= len(self.scoreList) - 1:
                            kFlag = 1
                            lFlag = 1
                        ### otherwise, only break out of inner loop:
                        #l = len(self.scoreList)
                        lFlag = 1
                    else:
                        l += 1
                        if l >= len(self.scoreList):
                            kFlag = 1
                            lFlag = 1
                        count += 1
            else:
                k += 1
                if k == (len(self.scoreList) - 1):
                    kFlag = 1
        '''

        ### evaluate pairs, if mismatched, update weight:

        ### if a record is in the top tau, it is supposed to be good

        self.traindp = (len(self.scoreList) * (self.tau / 100.0))
        self.trainnpr = 0

        #updates = 0
        i = 0
        iFlag = 0
        mistakes = 0

        #while i < ((len(self.scoreList) / 5) - 1):
        while iFlag == 0:
        #for i in range((len(self.scoreList) / 5) - 1):

            ### Only update if there is an irrelevant record
            ### ranked above a relevant record:

            if self.scoreList[i][1][1] == -1:
                #mistakes += 1 # comment out if counting mistakes by new method
                j = i + 1
                ### update if we have an irrelevant record scored
                ### higher than a relevant record:
                jFlag = 0
                while jFlag == 0:
                #for j in range (i + 1, (len(self.scoreList) / 5)):
                    if (self.scoreList[i][1][1] == self.scoreList[j][1][1]):
                        #mistakes += 1 # comment out if counting mistakes by new method
                        j += 1
                        if j >= len(self.scoreList) * self.tau / 100.0:
                            iFlag = 1
                            jFlag = 1
                    else:
                        self.eta = 1.0
                        self.updateWeight (self.scoreList[i], self.scoreList[j])
                        #updates += 1
                        #self.trainMistakes += mistakes # comment out if counting mistakes by new method
                        #mistakes = 0 # comment out if counting mistakes by new method
                        #i = j + 1 # update only using the first mistake
                        i += 1 # update every mistake
                        jFlag = 1
            else:
                i += 1
                self.trainnpr += 1
                if i >= len(self.scoreList) * self.tau / 100.0:
                    iFlag = 1


            '''
            else:
                ### Number of actual good in the top 20%
                self.trainnpr += 1
                    #else:
                        #self.trainnpr += 1
            '''
        #print 'updates: ', updates

    ### end pairwise function:

    #####
    ### newPairwise function:

    def newPairwise (self, sampleSizeK, thePath, gene, kmer, writeFile, logging):
        scores = []
        labels = []
        writeFile.write ('Train' + '\n')
        logging.info ('Train' + '\n')
        fileName = thePath + gene + '-' + str(kmer)

        # 1. sample K, k-mers:
        # 2. score k-mers using current weight
        # 3. sort them according to scores
        # 4. determine a decision threshold using desired precision
        # 5. using decision threshold, calculate:
        #    tp(fb) = records above threshold where label == 1
        #    fp(fb) = records above threshold where label == -1
        #    |Y+| (good*) = records in top tau where label == 1
        # 6. lower bound tp(fb) and upper bound fp(fb)
        #    as the zero-one loss:
        #    tp(fb) = sum records in top tau below threshold
        #       where label == 1
        #    fp(fb) = sum records above threshold where label == -1
        #    update weights for these records by:
        # update weight: w_t+1 = w_t + eta * (x_good - x_bad)

        ### 1. Sample k k-mers:

        ### create a list of size sampleSizeK of random line numbers
        ### from the file:
        sampleList = self.getRandomList(sampleSizeK)
        self.trainTotal = sampleSizeK
        self.trainMistakes = 0
        yPlus = 0

        ### 2. Score k-mers using current weights

        ### read the selected samples from the file into a list:
        ### and calculate their score:

        readFile = open (fileName, 'r')

        r = 0
        #self.scoreList = []
        scoreList = []

        for s in sampleList:
            t = s - r - 1
            ### iterate to the line before the one in s:
            for x in range(t):
                tempString = readFile.readline()
            ### read the sample:

            sampleX = processLine(readFile, 'bogus')
            if sampleX == 0:
                continue
            theScore = self.calcScore (sampleX)
            theLabel = sampleX[1]
            #self.scores.append(theScore)
            scores.append(theScore)

            ### labels are {-1, 1} but sklearn metrics uses {0,1}
            if theLabel == 1:
                #self.labels.append(theLabel)
                labels.append(theLabel)
                yPlus += 1
            else:
                #self.labels.append(0)
                labels.append(0)

            scoreTuple = (theScore, sampleX)
            #TODO **replace with insertSort
            #self.scoreList.append(scoreTuple)
            scoreList.append(scoreTuple)
            r = s

        readFile.close()

        ### 3. sort them according to scores

        ### sort the scores: (This could be replaced by insertSort)
        ### Put the highest scores at the top
        ### to rank them from highest to lowest
        #self.scoreList.sort(reverse=True)
        scoreList.sort(reverse=True)

        # 4. determine a decision threshold using desired precision

        #npScores = np.array(self.scores)
        #npLabels = np.array(self.labels)
        npScores = np.array(scores)
        npLabels = np.array(labels)
        positive_scores = npScores[npLabels == 1.0]
        #positive_scores = self.scores[self.labels == 1.0]

        #threshold = np.percentile(positive_scores, 100 - self.target_recall * 100)
        '''
        for tau in self.tauList:
            self.tau = tau
            print '\t\t\ttau: ', tau

            for target in self.targetList:
        '''
        writeFile.write(str(self.tau) + ',')
        logging.info(str(self.tau) + ',')
        #threshold = np.percentile(positive_scores, 100 - target * 100)
        self.threshold = np.percentile(positive_scores, self.target * 100)
        print 'threshold: ', self.threshold
        print 'target: ', self.target
        writeFile.write (str(self.threshold) + ',' + str(self.target) + ',')
        logging.info (str(self.threshold) + ',' + str(self.target) + ',')

        # 5. using decision threshold, calculate:
        #    tp(fb) = records above threshold where label == 1
        #    fp(fb) = records above threshold where label == -1
        #    |Y+| (good*) = records in top tau where label == 1

        tpfb = 0
        fpfb = 0
        #yPlus = 0
        tplfb = 0

        '''
        for j in range (len(self.scoreList)):
            if self.scoreList[j][1][1] == 1:
                yPlus += 1
        '''

        #for i in range (int(len(self.scoreList) * self.tau / 100)):
        for i in range (int(len(scoreList) * self.tau / 100)):
            ### if score is above threshold,
            ### if record is positive, increment true positives, tpfb
            ###   and increment |Y+| (yPlus)
            ### if record is negative, increment false positives, fpfb
            ### if score is below threshold and record is positive,
            ###   increment |Y+| (yPlus)

            '''
            if self.scoreList[i][0] >= self.threshold:
                if self.scoreList[i][1][1] == 1:
                    tpfb += 1
                    #yPlus += 1
                else:
                    fpfb += 1
            else:
                if self.scoreList[i][1][1] == 1:
                    #yPlus += 1
                    # if score is below threshold
                    # and label == 1, add to
                    # true positive lower bound of fb:
                    tplfb += 1
            '''

            if scoreList[i][0] >= self.threshold:
                if scoreList[i][1][1] == 1:
                    tpfb += 1
                else:
                    fpfb += 1
            else:
                if scoreList[i][1][1] == 1:
                    # if score is below threshold
                    # and label == 1, add to
                    # true positive lower bound of fb:
                    tplfb += 1

        print 'tpfb: ', tpfb
        print 'fpfb: ', fpfb
        print 'Y+: ', yPlus
        precision = float (tpfb) / (tpfb + fpfb)
        print '\t\tprecision: ', precision
        if yPlus > 0:
            recall = float(tpfb) / yPlus
        else:
            recall = 0
        print '\t\trecall: ', recall

        writeFile.write(str(tpfb) + ',')
        writeFile.write(str(fpfb) + ',')
        writeFile.write(str(yPlus) + ',')
        writeFile.write(str(precision) + ',')
        writeFile.write(str(recall) + '\n')
        logging.info(str(tpfb) + ',')
        logging.info(str(fpfb) + ',')
        logging.info(str(yPlus) + ',')
        logging.info(str(precision) + ',')
        logging.info(str(recall) + '\n')

        # 6. lower bound tp(fb) and upper bound fp(fb)
        #    as the zero-one loss:
        #    tp(fb) = sum records in top tau below threshold
        #       where label == 1
        #    fp(fb) = sum records above threshold where label == -1
        #    update weights for these records by:
        # update weight: w_t+1 = w_t + eta * (x_good - x_bad)

        ### evaluate pairs, if mismatched, update weight:

        ### if a record is in the top tau, it is supposed to be good

        #self.traindp = (len(self.scoreList) * (self.tau / 100.0))
        self.traindp = (len(scoreList) * (self.tau / 100.0))
        self.trainnpr = 0

        #updates = 0
        i = 0
        iFlag = 0
        mistakes = 0

        #while i < ((len(self.scoreList) / 5) - 1):
        while iFlag == 0:
        #for i in range((len(self.scoreList) / 5) - 1):

            ### Only update if there is an irrelevant record
            ### ranked above a relevant record:

            #if self.scoreList[i][1][1] == -1:
            if scoreList[i][1][1] == -1:
                #mistakes += 1 # comment out if counting mistakes by new method
                j = i + 1
                ### update if we have an irrelevant record scored
                ### higher than a relevant record:
                jFlag = 0
                while jFlag == 0:
                #for j in range (i + 1, (len(self.scoreList) / 5)):
                    #if (self.scoreList[i][1][1] == self.scoreList[j][1][1]):
                    if (scoreList[i][1][1] == scoreList[j][1][1]):
                        #mistakes += 1 # comment out if counting mistakes by new method
                        j += 1
                        #if j >= len(self.scoreList) * self.tau / 100.0:
                        if j >= len(scoreList) * self.tau / 100.0:
                            iFlag = 1
                            jFlag = 1
                    else:
                        self.eta = 1.0
                        #self.updateWeight (self.scoreList[i], self.scoreList[j])
                        self.updateWeight (scoreList[i], scoreList[j])
                        #updates += 1
                        #self.trainMistakes += mistakes # comment out if counting mistakes by new method
                        #mistakes = 0 # comment out if counting mistakes by new method
                        #i = j + 1 # update only using the first mistake
                        i += 1 # update every mistake
                        jFlag = 1
            else:
                i += 1
                self.trainnpr += 1
                #if i >= len(self.scoreList) * self.tau / 100.0:
                if i >= len(scoreList) * self.tau / 100.0:
                    iFlag = 1

    ### end newPairwise function:

    ### newPairwiseSGD function:

    def newPairwiseSGD (self, sampleSizeK, thePath, gene, kmer, writeFile, logging):
        writeFile.write ('Train' + '\n')
        logging.info ('Train' + '\n')
        fileName = thePath + gene + '-' + str(kmer)

        # 1. sample K, k-mers:
        # 2. score k-mers using current weight
        # 3. sort them according to scores
        # 4. determine a decision threshold using desired precision
        # 5. using decision threshold, calculate:
        #    tp(fb) = records above threshold where label == 1
        #    fp(fb) = records above threshold where label == -1
        #    |Y+| (good*) = records in top tau where label == 1
        # 6. lower bound tp(fb) and upper bound fp(fb)
        #    as the zero-one loss:
        #    tp(fb) = sum records in top tau below threshold
        #       where label == 1
        #    fp(fb) = sum records above threshold where label == -1
        #    update weights for these records by:
        # update weight: w_t+1 = w_t + eta * (x_good - x_bad)

        ### 1. Sample k k-mers:

        ### create a list of size sampleSizeK of random line numbers
        ### from the file:
        sampleList = self.getRandomList(sampleSizeK)
        self.trainTotal = sampleSizeK
        self.trainMistakes = 0

        ### 2. Score k-mers using current weights

        ### read the selected samples from the file into a list:
        ### and calculate their score:

        readFile = open (fileName, 'r')

        r = 0
        self.scoreList = []

        for s in sampleList:
            t = s - r - 1
            ### iterate to the line before the one in s:
            for x in range(t):
                tempString = readFile.readline()
            ### read the sample:
            sampleX = processLine(readFile, 'bogus')
            if sampleX == 0:
                continue
            theScore = self.calcScore (sampleX)
            theLabel = sampleX[1]
            self.scores.append(theScore)

            ### labels are {-1, 1} but sklearn metrics uses {0,1}
            if theLabel == 1:
                self.labels.append(theLabel)
            else:
                self.labels.append(0)

            scoreTuple = (theScore, sampleX)
            #TODO **replace with insertSort
            self.scoreList.append(scoreTuple)
            r = s

        readFile.close()

        ### 3. sort them according to scores

        ### sort the scores: (This could be replaced by insertSort)
        ### Put the highest scores at the top
        ### to rank them from highest to lowest
        self.scoreList.sort(reverse=True)

        # 4. determine a decision threshold using desired precision

        npScores = np.array(self.scores)
        npLabels = np.array(self.labels)
        positive_scores = npScores[npLabels == 1.0]
        #positive_scores = self.scores[self.labels == 1.0]

        #threshold = np.percentile(positive_scores, 100 - self.target_recall * 100)

        for tau in self.tauList:
            self.tau = tau
            print '\t\t\ttau: ', tau

            for target in self.targetList:
                writeFile.write(str(tau) + ',')
                logging.info(str(tau) + ',')
                #threshold = np.percentile(positive_scores, 100 - target * 100)
                self.threshold = np.percentile(positive_scores, target * 100)
                print 'threshold: ', self.threshold
                print 'target: ', target
                writeFile.write (str(self.threshold) + ',' + str(target) + ',')
                logging.info (str(self.threshold) + ',' + str(target) + ',')

                # 5. using decision threshold, calculate:
                #    tp(fb) = records above threshold where label == 1
                #    fp(fb) = records above threshold where label == -1
                #    |Y+| (good*) = records in top tau where label == 1

                tpfb = 0
                fpfb = 0
                yPlus = 0
                tplfb = 0

                for i in range (int(len(self.scoreList) * self.tau / 100)):
                    ### if score is above threshold,
                    ### if record is positive, increment true positives, tpfb
                    ###   and increment |Y+| (yPlus)
                    ### if record is negative, increment false positives, fpfb
                    ### if score is below threshold and record is positive,
                    ###   increment |Y+| (yPlus)

                    if self.scoreList[i][0] >= self.threshold:
                        if self.scoreList[i][1][1] == 1:
                            tpfb += 1
                            yPlus += 1
                        else:
                            fpfb += 1
                    else:
                        if self.scoreList[i][1][1] == 1:
                            yPlus += 1
                            # if score is below threshold
                            # and label == 1, add to
                            # true positive lower bound of fb:
                            tplfb += 1

                print 'tpfb: ', tpfb
                print 'fpfb: ', fpfb
                print 'Y+: ', yPlus
                precision = float (tpfb) / (tpfb +fpfb)
                print '\t\tprecision: ', precision
                if yPlus > 0:
                    recall = float(tpfb) / yPlus
                else:
                    recall = 0
                print '\t\trecall: ', recall

                writeFile.write(str(tpfb) + ',')
                writeFile.write(str(fpfb) + ',')
                writeFile.write(str(yPlus) + ',')
                writeFile.write(str(precision) + ',')
                writeFile.write(str(recall) + '\n')
                logging.info(str(tpfb) + ',')
                logging.info(str(fpfb) + ',')
                logging.info(str(yPlus) + ',')
                logging.info(str(precision) + ',')
                logging.info(str(recall) + '\n')

            # 6. lower bound tp(fb) and upper bound fp(fb)
            #    as the zero-one loss:
            #    tp(fb) = sum records in top tau below threshold
            #       where label == 1
            #    fp(fb) = sum records above threshold where label == -1
            #    update weights for these records by:
            # update weight: w_t+1 = w_t + eta * (x_good - x_bad)

            for i in range (int(len(self.scoreList) * self.tau / 100)):
                '''
                if self.scoreList[i][0] >= self.threshold:
                    if self.scoreList[i][1][1] == -1:
                        self.SGDUpdate (self.scoreList[i])

                else:
                    if self.scoreList[i][1][1] == 1:
                        #for upSample in range(4):
                        self.SGDUpdate (self.scoreList[i])
                '''
                if self.scoreList[i][0] < self.threshold:
                    if self.scoreList[i][1][1] == 1:
                        self.SGDUpdate (self.scoreList[i])


    ### end newPairwiseSGD function:

    def SGDUpdate (self, i):
        x = i[1][0]
        yStar = i[1][1]

        for index in x:
            #print 'index: ', index[0]
            #print 'value: ', index[1]
            self.w[index[0]] += self.eta * yStar * index[1]

    ### end SGDUpdate


    #####

    ### newTest function:

    def newTest (self, sampleSizeK, thePath, gene, kmer, writeFile, logging):
        scores = []
        labels = []
        writeFile.write ('Test:' + '\n')
        logging.info ('Test:' + '\n')
        fileName = thePath + gene + '-' + str(kmer)

        # 2. score k-mers using current weight
        # 3. sort them according to scores
        # 4. determine a decision threshold using desired precision
        # 5. using decision threshold, calculate:
        #    tp(fb) = records above threshold where label == 1
        #    fp(fb) = records above threshold where label == -1
        #    |Y+| (good*) = records in top tau where label == 1
        # 6. lower bound tp(fb) and upper bound fp(fb)
        #    as the zero-one loss:
        #    tp(fb) = sum records in top tau below threshold
        #       where label == 1
        #    fp(fb) = sum records above threshold where label == -1
        #    update weights for these records by:
        # update weight: w_t+1 = w_t + eta * (x_good - x_bad)

        self.testTotal = 0
        self.testMistakes = 0

        ### 2. Score k-mers using current weights

        ### read the selected samples from the file into a list:
        ### and calculate their score:

        readFile = open (fileName, 'r')

        r = 0
        #self.scoreList = []
        scoreList = []

        lineCount = 0
        yPlus = 0
        #tempString = readFile.readline()
        sampleX = processLine(readFile, 'bogus')

        #while tempString:
        while sampleX != 0:
            lineCount += 1
            self.testTotal += 1
            ### read the sample:
            #sampleX = processLine(readFile, 'bogus')
            '''
            if sampleX == 0:
                continue
            '''
            theScore = self.calcScore (sampleX)
            theLabel = sampleX[1]
            #self.scores.append(theScore)
            scores.append(theScore)

            ### labels are {-1, 1} but sklearn metrics uses {0,1}
            if theLabel == 1:
                #self.labels.append(theLabel)
                labels.append(theLabel)
                yPlus += 1
            else:
                #self.labels.append(0)
                labels.append(0)

            scoreTuple = (theScore, sampleX)
            #TODO **replace with insertSort
            #self.scoreList.append(scoreTuple)
            scoreList.append(scoreTuple)

            #tempString = readFile.readline()
            sampleX = processLine(readFile, 'bogus')

        readFile.close()
        print '\t\t***\t*** lineCount: ', lineCount, ' ***\t***'

        ### 3. sort them according to scores

        ### sort the scores: (This could be replaced by insertSort)
        ### Put the highest scores at the top
        ### to rank them from highest to lowest
        #self.scoreList.sort(reverse=True)
        scoreList.sort(reverse=True)

        # 4. determine a decision threshold using desired precision

        #npScores = np.array(self.scores)
        #npLabels = np.array(self.labels)
        npScores = np.array(scores)
        npLabels = np.array(labels)
        positive_scores = npScores[npLabels == 1.0]
        #positive_scores = self.scores[self.labels == 1.0]

        #threshold = np.percentile(positive_scores, 100 - self.target_recall * 100)

        '''
        for tau in self.tauList:
            self.tau = tau
            print '\t\t\ttau: ', tau

            for target in self.targetList:
        '''
        writeFile.write(str(self.tau) + ',')
        logging.info(str(self.tau) + ',')
        #threshold = np.percentile(positive_scores, 100 - target * 100)
        self.threshold = np.percentile(positive_scores, self.target * 100)
        print 'threshold: ', self.threshold
        print 'target: ', self.target
        writeFile.write (str(self.threshold) + ',' + str(self.target) + ',')
        logging.info (str(self.threshold) + ',' + str(self.target) + ',')

        # 5. using decision threshold, calculate:
        #    tp(fb) = records above threshold where label == 1
        #    fp(fb) = records above threshold where label == -1
        #    |Y+| (good*) = records in top tau where label == 1

        tpfb = 0
        fpfb = 0
        #yPlus = 0
        tplfb = 0

        '''
        ### yPlus is the count of all records with a label == 1:
        for i in range (len(self.scoreList)):
            if self.scoreList[i][1][1] == 1:
                yPlus += 1
        '''

        #for i in range (int(len(self.scoreList) * self.tau / 100)):
        for i in range (int(len(scoreList) * self.tau / 100)):
            ### if score is above threshold,
            ### if record is positive, increment true positives, tpfb
            ###   and increment |Y+| (yPlus)
            ### if record is negative, increment false positives, fpfb
            ### if score is below threshold and record is positive,
            ###   increment |Y+| (yPlus)

            '''
            if self.scoreList[i][0] >= self.threshold:
                if self.scoreList[i][1][1] == 1:
                    tpfb += 1
                    #yPlus += 1
                else:
                    fpfb += 1
            else:
                if self.scoreList[i][1][1] == 1:
                    #yPlus += 1
                    # if score is below threshold
                    # and label == 1, add to
                    # true positive lower bound of fb:
                    tplfb += 1
            '''

            if scoreList[i][0] >= self.threshold:
                if scoreList[i][1][1] == 1:
                    tpfb += 1
                else:
                    fpfb += 1
            else:
                if scoreList[i][1][1] == 1:
                    # if score is below threshold
                    # and label == 1, add to
                    # true positive lower bound of fb:
                    tplfb += 1
        
        print 'tpfb: ', tpfb
        print 'fpfb: ', fpfb
        print 'Y+: ', yPlus
        precision = float (tpfb) / (tpfb +fpfb)
        print '\t\tprecision: ', precision
        recall = float(tpfb) / yPlus
        print '\t\trecall: ', recall

        writeFile.write(str(tpfb) + ',')
        writeFile.write(str(fpfb) + ',')
        writeFile.write(str(yPlus) + ',')
        writeFile.write(str(precision) + ',')
        writeFile.write(str(recall) + '\n')
        logging.info(str(tpfb) + ',')
        logging.info(str(fpfb) + ',')
        logging.info(str(yPlus) + ',')
        logging.info(str(precision) + ',')
        logging.info(str(recall) + '\n')



    ### end newTest function:

    #def testWeight (self, trainDataSet, testDataSet, outFile, t):
    def testWeight (self, testDataSet, t):
        self.testTotal += len(testDataSet)

        xit = [] # x-sub-i-sub-t, the training vector for the current iterations
        yit = 0 # y-sub-i-sub-t, the label for the current training vector (yStar)
        yHat = 0.0

        ### evaluate all test samples:
        #self.testTotal = 0
        #self.testMistakes = 0

        for i in range(len(testDataSet)):
            xit = testDataSet[i][0]
            yit = testDataSet[i][1]
            if yit == 1:
                self.testGood += 1 ### increment count of actual good

            ### make the prediction:
            #yHat = yit * self.dotProd(xit) # Method A
            #yHat = yit * (self.dotProd(xit) + self.b) # Method B
            #yHat = sign(yit * (self.dotProd(xit) + self.b)) # Method C
            #yHat = sign(yit * self.dotProd(xit)) # Method D
            yHat = sign(self.dotProd(xit)) # Correction

            ### if the value is actually good,
            ### and we predicted good, increment npr
            ### which is the # lines predicted as good that are actually good
            if yHat > 0: # if we predicted good:
                self.testdp += 1 # +1 prediction count
                if yit == 1: # if it is actually good:
                    self.testnpr += 1 # correct +1 prediction
                else:
                    self.testMistakes += 1
            else: ### yHat <= 0: # if we predicted bad:
                if yit == 1: # if it is actually good:
                    self.testMistakes += 1

        print 'iteration: ', t
        print 'testMistakes: ', self.testMistakes
        #print 'testTotal: ', self.testTotal

    ### end testWeight function

    ### testWeightPair

    def testWeightPair (self, testDataSet, t):
        self.testTotal += len(testDataSet) * self.tau / 100.0
        #self.testMistakes = 0
        #self.testGood = 0
        testScoreList = []

        for feature in testDataSet:
            theScore = self.calcScore(feature)
            testScoreList.append([theScore, feature])
        testScoreList.sort(reverse=True)

        #print testScoreList[0]
        #print testScoreList[len(testScoreList) / 5 - 1]

        '''
        ### count total good records:
        for feature in range(int(len(testDataSet) / self.tau * 100 )):
            if testDataSet[feature][1] == 1:
                self.testGood += 1
        '''
        #for k in range(int (len(testScoreList) * self.tau / 100.0)):
        for k in range(len(testScoreList)):
            if self.scoreList[k][1][1] == 1:
                self.testgoodStar += 1

        for x in range (int (len(testScoreList) * self.tau / 100.0)):
            if testScoreList[x][1][1] == -1:
                self.testMistakes += 1

        self.testdp = len(testScoreList) * self.tau / 100.0

        '''
        ### a mistake is an irrelevant record scored higher than a relevant record:
        ### traverse the list and count irrelevant records until a relevant record
        ### is located then add that to the mistakes sum:

        k = 0
        flag = 0
        goodIndex = []
        while flag == 0:
            count = 1
            if testScoreList[k][1][1] == -1:
                l = k + 1
                while (l < len(testScoreList) / 5):
                    if (testScoreList[l][1][1] == 1):
                        goodIndex.append(l)
                        self.testMistakes += count
                        k = l
                        if k == len(testScoreList) / 5 - 1:
                            flag = 1
                        l = len(testScoreList)
                    else:
                        l += 1
                        if l == len(testScoreList) / 5:
                            flag = 1
                        count += 1
            else:
                k += 1
                if k == (len(testScoreList) / 5 - 1):
                    flag = 1

            '''


        print 'iteration: ', t
        print 'testMistakes: ', self.testMistakes
        #print goodIndex
    ### end testWeightPair function

    ### function dotProd(), dot product by index
    def dotProd(self, xArray):
        result = 0.0

        for element in xArray:
            if self.w.get(element[0], '--') == '--':
                #self.w[element[0]] = 0
                continue
            else:
                result += self.w[element[0]] * element[1]

        return result

    ### end dotProd() function

    ### saveWeight function, save the weight vector to file:

    def saveWeight (self):
        writeFile = open ('weight.txt', 'w')
        weightList = self.w.keys()
        weightList.sort()
        for index in weightList:
            writeFile.write ( str(index) + ', ' + str(self.w.get(index)) + '\n')
        writeFile.close()

    ### end saveWeight function

    ### get a list of random numbers of size sampleSizeK

    def getRandomList (self, sampleSizeK):
        randomList = []
        for s in range (sampleSizeK):
            randomList.append( random.randrange(0, self.fileLength) )
        randomList.sort()
        return randomList

    ### end getRandomList function

    ### calcScore function, calculate the score of a sample k-mer:

    def calcScore (self, sampleX):
        xit = sampleX[0]
        yit = sampleX[1]
        theScore = 0.0

        ### The score is the dot-product of the weight and the
        ### feature vector:

        for index in xit:
            if self.w.get(index[0], '--') == '--':
                self.w[index[0]] = 0.0
            theScore += self.w.get(index[0]) * index[1] * yit

        return theScore

    ### end calcScore function

    ### testcalcScore function, calculate the score of a sample k-mer:
    ### What if we only use weights that match the sign of our label?

    def testcalcScore (self, sampleX):
        xit = sampleX[0]
        yit = sampleX[1]
        theScore = 0.0

        ### The score is the dot-product of the weight and the
        ### feature vector:

        for index in xit:
            if self.w.get(index[0], '--') == '--':
                self.w[index[0]] = 0.0
            if sign(self.w.get(index[0])) == sign(yit):
                theScore += self.w.get(index[0]) * index[1] * yit

        return theScore

    ### end testcalcScore function

    ### update the weight:


    def updateWeight (self, i, j):
        ### i is the irrelevant record,
        ### subtract its weights:
        ### j is the relevant record,
        ### add its weights:

        '''
        if self.underSample == 3:
            for index in i[1][0]:
                self.w[index[0]] -= self.eta * index[1]
                self.underSample = 0
        else:
            self.underSample += 1



        for index in i[1][0]:
            self.w[index[0]] -= self.eta * index[1]
        '''

        for index in j[1][0]:
            self.w[index[0]] += self.eta * index[1]
            #self.w[index[0]] += self.eta * index[1]
        '''
        ### This is an expensive evaluation that looks at each index from each feature set
        ### but the above achieves the same result with less complexity
        ### (unless we want a different learning rate for each classifier)

        ### initialize indexes for each feature:
        ii = 0
        jj = 0

        ### Check each index in each feature set,
        ### If they match, use the difference to update,
        ### otherwise, use the value of the unmatched index

        iFlag = 0
        jFlag = 0
        while iFlag * jFlag != 1:

            if iFlag == 0:
                if jFlag == 0:
                    ### if the current index of i is lower than the current index of j:
                    if (i[1][0][ii][0] < j[1][0][jj][0]):
                        self.w[i[1][0][ii][0]] = self.w[i[1][0][ii][0]] - self.eta * i[1][0][ii][1]
                        #self.w[i[1][0][ii][0]] = self.w[i[1][0][ii][0]] - self.eta * i[1][0][ii][1] * i[1][1]
                        ii += 1
                        if ii == len(i[1][0]):
                            iFlag = 1
                    ### if the current index of i matches the current index of j:
                    elif (i[1][0][ii][0] == j[1][0][jj][0]):
                        #self.w[i[1][0][ii][0]] = self.w[i[1][0][ii][0]] + self.eta * ((j[1][0][jj][0] - i[1][0][ii][0] ))
                        #self.w[i[1][0][ii][0]] = self.w[i[1][0][ii][0]] + self.eta * ((i[1][0][ii][0] * i[1][1]) + (j[1][0][jj][0] * j[1][1]))
                        ii += 1
                        jj += 1
                        if ii == len(i[1][0]):
                            iFlag = 1
                        if jj == len(j[1][0]):
                            jFlag = 1
                    ### if the current index of i is higher than the current index of j:
                    elif (i[1][0][ii][0] > j[1][0][jj][0]):
                        self.w[j[1][0][jj][0]] = self.w[j[1][0][jj][0]] + self.eta * j[1][0][jj][1]
                        #self.w[j[1][0][jj][0]] = self.w[j[1][0][jj][0]] + self.eta * j[1][0][jj][1] * j[1][1]
                        jj += 1
                        if jj == len(j[1][0]):
                            jFlag = 1
                ### jFlag == 1, there are no more indexes in j:
                else:
                    self.w[i[1][0][ii][0]] = self.w[i[1][0][ii][0]] - self.eta * i[1][0][ii][1]
                    #self.w[i[1][0][ii][0]] = self.w[i[1][0][ii][0]] + self.eta * i[1][0][ii][1] * i[1][1]
                    ii += 1
                    if ii == len(i[1][0]):
                        iFlag = 1

            ### iFlag == 1, no more indexes in i, iterate through remainder of j:
            else:
                if (jFlag == 0):
                    self.w[j[1][0][jj][0]] = self.w[j[1][0][jj][0]] + self.eta  * j[1][0][jj][1]
                    #self.w[j[1][0][jj][0]] = self.w[j[1][0][jj][0]] + self.eta * j[1][0][jj][1] * j[1][1]
                    jj += 1
                    if jj == len(j[1][0]):
                        jFlag = 1
        '''

        '''
            if iFlag == 0:
                if jFlag == 0:
                    ### if the current index of i is lower than the current index of j:
                    if (i[1][0][ii][0] < j[1][0][jj][0]):
                        if self.w[i[1][0][ii][0]] != 0:
                            if sign(self.w[i[1][0][ii][0]]) == sign(i[1][1]):
                                self.w[i[1][0][ii][0]] = self.w[i[1][0][ii][0]] - self.eta * i[1][0][ii][1]
                            else:
                                self.w[i[1][0][ii][0]] = 0
                        else:
                            self.w[i[1][0][ii][0]] = self.w[i[1][0][ii][0]] - self.eta * i[1][0][ii][1]
                        #self.w[i[1][0][ii][0]] = self.w[i[1][0][ii][0]] - self.eta * i[1][0][ii][1] * i[1][1]
                        ii += 1
                        if ii == len(i[1][0]):
                            iFlag = 1
                    ### if the current index of i matches the current index of j:
                    elif (i[1][0][ii][0] == j[1][0][jj][0]):
                        #self.w[i[1][0][ii][0]] = self.w[i[1][0][ii][0]] + self.eta * ((j[1][0][jj][0] - i[1][0][ii][0] ))
                        #self.w[i[1][0][ii][0]] = self.w[i[1][0][ii][0]] + self.eta * ((i[1][0][ii][0] * i[1][1]) + (j[1][0][jj][0] * j[1][1]))
                        ii += 1
                        jj += 1
                        if ii == len(i[1][0]):
                            iFlag = 1
                        if jj == len(j[1][0]):
                            jFlag = 1
                    ### if the current index of i is higher than the current index of j:
                    elif (i[1][0][ii][0] > j[1][0][jj][0]):
                        if self.w[j[1][0][jj][0]] != 0:
                            if sign(self.w[j[1][0][jj][0]]) == sign(j[1][1]):
                                self.w[j[1][0][jj][0]] = self.w[j[1][0][jj][0]] + self.eta * j[1][0][jj][1]
                            else:
                                self.w[j[1][0][jj][0]] = 0
                        else:
                            self.w[j[1][0][jj][0]] = self.w[j[1][0][jj][0]] + self.eta * j[1][0][jj][1]
                        #self.w[j[1][0][jj][0]] = self.w[j[1][0][jj][0]] + self.eta * j[1][0][jj][1] * j[1][1]
                        jj += 1
                        if jj == len(j[1][0]):
                            jFlag = 1
                ### jFlag == 1, there are no more indexes in j:
                else:
                    if self.w[i[1][0][ii][0]] != 0:
                        if sign(self.w[i[1][0][ii][0]]) == sign(i[1][1]):
                            self.w[i[1][0][ii][0]] = self.w[i[1][0][ii][0]] - self.eta * i[1][0][ii][1]
                        else:
                            self.w[i[1][0][ii][0]] = 0
                    else:
                        self.w[i[1][0][ii][0]] = self.w[i[1][0][ii][0]] - self.eta * i[1][0][ii][1]
                    #self.w[i[1][0][ii][0]] = self.w[i[1][0][ii][0]] + self.eta * i[1][0][ii][1] * i[1][1]
                    ii += 1
                    if ii == len(i[1][0]):
                        iFlag = 1

            ### iFlag == 1, no more indexes in i, iterate through remainder of j:
            else:
                if (jFlag == 0):
                    if self.w[j[1][0][jj][0]] != 0:
                        if sign(self.w[j[1][0][jj][0]]) == sign(j[1][1]):
                            self.w[j[1][0][jj][0]] = self.w[j[1][0][jj][0]] + self.eta  * j[1][0][jj][1]
                        else:
                            self.w[j[1][0][jj][0]] = 0
                    else:
                        self.w[j[1][0][jj][0]] = self.w[j[1][0][jj][0]] + self.eta  * j[1][0][jj][1]
                    #self.w[j[1][0][jj][0]] = self.w[j[1][0][jj][0]] + self.eta * j[1][0][jj][1] * j[1][1]
                    jj += 1
                    if jj == len(j[1][0]):
                        jFlag = 1
        '''

    ### end updateWeight function

### End of file
