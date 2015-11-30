import math
import numpy as np

class CutOff:
    def __init__(self, posScores, negScores):
        # Total Scores
        self.totalPosScores = len(posScores)
        self.totalNegScores = len(negScores)

        # Sort the Scores
        self.posScoresList = posScores
        self.negScoresList = negScores

        self.falsePositiveRateList = None
        self.truePositiveRateList = None

        self.value = None
        self.sensitivityPercent = None
        self.specificityPercent = None
        self.mcc = None

    def calculatePositive(self, scoreList, probableCutOff):
        count = 0
        for i in scoreList:
            if i > probableCutOff:
                count +=1
        return count

    def calcSensitivity(self):
        truePositive = self.calculatePositive(self.posScoresList, self.value)
        self.sensitivityPercent = float(truePositive) / self.totalPosScores * 100

    def calcSpecificity(self):
        trueNegative=0
        falsePositive=0

        for i in self.negScoresList:
            if i <= self.value:
                trueNegative+=1

        for i in self.posScoresList:
            if i <= self.value:
                falsePositive+=1

        self.specificityPercent = float(trueNegative) / (falsePositive + trueNegative) * 100

    def getCutOff(self, probableCutOffStart, probableCutOffEnd):
        print "calculating cut off .... \t\t\t " ,
        self.value = 0
        currentMcc = 0
        
        if min(self.posScoresList) >= max(self.negScoresList):
            return min(self.posScoresList)

        # for i in np.arange(probableCutOffStart, probableCutOffEnd, 0.01):
        for i in np.arange(3, 4, 0.01):
            probableCutOff = i

            truePositive = self.calculatePositive(self.posScoresList, probableCutOff)
            falsePositive = self.calculatePositive(self.negScoresList, probableCutOff)

            trueNegative = self.totalNegScores - falsePositive
            falseNegative = self.totalPosScores - truePositive

            newMcc = float(truePositive * trueNegative - falsePositive * falseNegative) / \
                     (math.sqrt((truePositive + falsePositive) * (trueNegative + falsePositive) *
                                (truePositive + falseNegative) * (trueNegative + falseNegative))
                      )

            if currentMcc < newMcc:
                self.value = probableCutOff
                currentMcc = newMcc
                self.mcc = currentMcc

        self.calcSensitivity()
        self.calcSpecificity()
        print "done!"