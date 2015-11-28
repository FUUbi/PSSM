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

        self.sensitivity = None
        self.specificity = None

    def calculateCutOff(self, start, stop):
        print "calculating cut off .... \t\t\t " ,
        cutOffValue = 0
        currentMcc = 0

        for i in np.arange(start, stop, 0.1):
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
                cutOffValue = probableCutOff
                currentMcc = newMcc
        print "done!"
        return cutOffValue

    def calculatePositive(self, scoreList, probableCutOff):
        anzahl = 0
        for i in scoreList:
            if i > probableCutOff:
                anzahl +=1
        return anzahl

    def getCutOff(self, probableCutOffStart, probableCutOffEnd):
        cutOff = self.calculateCutOff(probableCutOffStart, probableCutOffEnd)
        return cutOff