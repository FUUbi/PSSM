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
    def probCutOff(self, start, stop, step):
        probCutOffList = list()
        i = start
        while i <= stop:
            i = round(i, 1)
            probCutOffList.append(i)
            i+= step
        return probCutOffList


    ## The TruePossitives "TP" for a Probable Cut Off Value,ScoreList "posScores"
    def calculateTruePositive(self, pco, posScores):
        truePositiveList = list()
        pco = pco
        posScores = posScores
        value = 0
        for i in range(len(pco)):
            for x in range(len(posScores)):
                if posScores[x] > pco[i]:
                    value +=1
            truePositiveList.append(value)
            value = 0
        return truePositiveList


    ## The TrueNegatives "TN" for a Probable Cut Off Value, ScoreList "negScores"
    def calculateTrueNegative(self, probCutOff, negScores):
        trueNegativeList = list()
        probCutOff = probCutOff
        negScores = negScores
        value = 0
        for i in range(len(probCutOff)):
            for x in range(len(negScores)):
                if negScores[x] < probCutOff[i]:
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


    def setCutOff(self, probCutOffList):
        tmin = math.sqrt(self.falsePositiveRateList[len(self.falsePositiveRateList) - 1] ** 2 + (self.truePositiveRateList[len(self.falsePositiveRateList) - 1] - 1) ** 2)
        for i in range(len(probCutOffList)-1, -1, -1):
            t = math.sqrt(self.falsePositiveRateList[i] ** 2 + (self.truePositiveRateList[i] - 1) ** 2)
            if t < tmin:
                self.cutOffValue = probCutOffList[i]
                self.sensitivity = self.falsePositiveRateList[i]
                self.specificity = self.truePositiveRateList[i]


    def rocCurve(self, probCutOffStart, probCutOffEnd, probCutOffStep):
        print("generating the Data for the ROC Curve")
        probCutOffStart = probCutOffStart
        probCutOffEnd = probCutOffEnd
        probCutOffStep = probCutOffStep

        probCutOffList = self.probCutOff(probCutOffStart, probCutOffEnd, probCutOffStep)
        tpList = self.calculateTruePositive(probCutOffList, self.posScoresSorted)
        tnList = self.calculateTrueNegative(probCutOffList, self.negScoresSorted)
        tnRateList = self.calculateTrueNegativeRate(tnList, self.totalNegScores)
        self.falsePositiveRateList = self.calculateFalsePositiveRate(tnRateList)
        self.truePositiveRateList = self.truePositiveRate(tpList, self.totalPosScores)

        self.setCutOff(probCutOffList)

