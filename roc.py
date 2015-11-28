import math


class Roc:
    def __init__(self, posScores, negScores):
        # Total Scores
        self.totalPosScores = len(posScores)
        self.totalNegScores = len(negScores)

        # Sort the Scores
        self.posScoresSorted = sorted(posScores, reverse=True)
        self.negScoresSorted = sorted(negScores, reverse=True)

        self.falsePositiveRateList = None
        self.truePositiveRateList = None

        self.cutOffValue = None
        self.sensitivity = None
        self.specificity = None


    ## Generate a List of Probable Cut Off Values
    ## start = start Value, stop = stop Value, step = steps between the Values
    def probableCutOff(self, start, stop, step):
        probableCutOffList = list()
        i = start
        while i <= stop:
            i = round(i, 1)
            probableCutOffList.append(i)
            i+= step
        return probableCutOffList


    ## The TruePossitives "TP" for a Probable Cut Off Value,ScoreList "posScores"
    def calculateTruePositive(self, probableCutOff, posScores):
        truePositiveList = list()
        probableCutOff = probableCutOff
        posScores = posScores
        value = 0
        for i in range(len(probableCutOff)):
            for x in range(len(posScores)):
                if posScores[x] > probableCutOff[i]:
                    value +=1
            truePositiveList.append(value)
            value = 0
        return truePositiveList


    ## The TrueNegatives "TN" for a Probable Cut Off Value, ScoreList "negScores"
    def calculateTrueNegative(self, probableCutOff, negScores):
        trueNegativeList = list()
        probableCutOff = probableCutOff
        negScores = negScores
        value = 0
        for i in range(len(probableCutOff)):
            for x in range(len(negScores)):
                if negScores[x] < probableCutOff[i]:
                    value +=1
            trueNegativeList.append(value)
            value = 0
        return trueNegativeList

    ## The TN-Rate "Specifity" i = length of the tnList, tnList = list of TN Values, totalNeg = total negative values
    def calculateTrueNegativeRate(self, trueNegativeList, totalNeg):
        trueNegativRateList = list()
        trueNegativeList = trueNegativeList
        totalNeg = totalNeg
        for i in range(len(trueNegativeList)):
            value = float(trueNegativeList[i]) / totalNeg
            trueNegativRateList.append(value)
        return trueNegativRateList

    ## The FP-Rate "1-TN Rate"
    def calculateFalsePositiveRate(self, trueNegativeRateList):
        falsePositiveRateList = list()
        trueNegativeRateList = trueNegativeRateList
        for i in range(len(trueNegativeRateList)):
            value = float(1 - trueNegativeRateList[i])
            falsePositiveRateList.append(value)
        return falsePositiveRateList

    ## The TP-Rate "Sensivity" tpList = list of TP Values, totalPos = total positive values
    def truePositiveRate(self, truePositiveList, totalPos):
        truePositiveRateList = list()
        truePositiveList = truePositiveList
        totalPos = totalPos
        for i in range(len(truePositiveList)):
            value = float(truePositiveList[i]) / totalPos
            truePositiveRateList.append(value)
        return truePositiveRateList


    def setCutOff(self, probableCutOffList):
        tmin = math.sqrt(self.falsePositiveRateList[len(self.falsePositiveRateList) - 1] ** 2 + (self.truePositiveRateList[len(self.falsePositiveRateList) - 1] - 1) ** 2)
        for i in range(len(probableCutOffList)-1, -1, -1):
            t = math.sqrt(self.falsePositiveRateList[i] ** 2 + (self.truePositiveRateList[i] - 1) ** 2)
            if t < tmin:
                self.cutOffValue = probableCutOffList[i]
                self.sensitivity = self.falsePositiveRateList[i]
                self.specificity = self.truePositiveRateList[i]


    def rocCurve(self, probableCutOffStart, probableCutOffEnd, probableCutOffStep):
        print("generating the Data for the ROC Curve")
        probableCutOffStart = probableCutOffStart
        probableCutOffEnd = probableCutOffEnd
        probableCutOffStep = probableCutOffStep

        probableCutOffList = self.probableCutOff(probableCutOffStart, probableCutOffEnd, probableCutOffStep)
        truePositiveList = self.calculateTruePositive(probableCutOffList, self.posScoresSorted)
        trueNegativeList = self.calculateTrueNegative(probableCutOffList, self.negScoresSorted)
        trueNegativeRateList = self.calculateTrueNegativeRate(trueNegativeList, self.totalNegScores)
        self.falsePositiveRateList = self.calculateFalsePositiveRate(trueNegativeRateList)
        self.truePositiveRateList = self.truePositiveRate(truePositiveList, self.totalPosScores)

        self.setCutOff(probableCutOffList)

