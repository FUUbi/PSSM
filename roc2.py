import math
import numpy as np


class CutOff:
    def __init__(self, cutOffValue, truePositive, falsePositive,trueNegative, falseNegative ):
        self.value = cutOffValue
        self.sensitivity = float(truePositive / (truePositive + falseNegative))
        self.specificity = float(trueNegative / (trueNegative + falsePositive))
        print self.sensitivity
        self.matthewscc = float(truePositive * trueNegative - falsePositive * falseNegative) / \
                                    (math.sqrt((truePositive + falsePositive) * (trueNegative + falsePositive) *
                                               (truePositive + falseNegative) * (trueNegative + falseNegative))
                                     )


class Roc2:
    def __init__(self, posScores, negScores):
        # Total Scores
        self.totalPosScores = len(posScores)
        self.totalNegScores = len(negScores)

        # Sort the Scores
        self.posScoresList = posScores
        self.negScoresList = negScores

        self.cutOffList = list()

        self.cutOff = None



    def calculatePositive(self, scoreList, probableCutOff):
        anzahl = 0
        for i in scoreList:
            if i > probableCutOff:
                anzahl +=1
        return anzahl

    def setCutOff(self, probableCutOffStart, probableCutOffEnd):
        print "calculating cut off .... \t\t\t " ,
        for probableCutOff in np.arange(probableCutOffStart, probableCutOffEnd, 0.1):

            truePositive = self.calculatePositive(self.posScoresList, probableCutOff)
            falsePositive = self.calculatePositive(self.negScoresList, probableCutOff)

            trueNegative = self.totalNegScores - falsePositive
            falseNegative = self.totalPosScores - truePositive
            newCutOff = CutOff(probableCutOff, truePositive, falsePositive, trueNegative, falseNegative)

            if(self.cutOff == None):
                self.cutOff = newCutOff

            if newCutOff.matthewscc > self.cutOff.matthewscc:
                self.cutOff = newCutOff

            self.cutOffList.append(newCutOff)

        print "done!"

    # def getSensitivity(self):
    #     sensitivities = []
    #     for cutOff in self.cutOffList:
    #         sensitivities.append(cutOff.sensitivity)
    #         print cutOff.matthewscc
    #     return sensitivities
    #
    # def getSpecificity(self):
    #     specificities = []
    #     for cutOff in self.cutOffList:
    #         print cutOff.sensitivity
    #         specificities.append(cutOff.sensitivity)
    #     return specificities