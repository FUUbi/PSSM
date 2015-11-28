import math


class Roc:
    def __init__(self, posScores, negScores):
        # Total Scores
        self.totalPosScores = len(posScores)
        self.totalNegScores = len(negScores)

        # Sort the Scores
        self.posScoresSorted = sorted(posScores, reverse=True)
        self.negScoresSorted = sorted(negScores, reverse=True)

        self.falsePositivRateList = None
        self.truePositivRateList = None

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
    def calculateTruePositives(self, pco, posScores):
        truePositivList = list()
        pco = pco
        posScores = posScores
        value = 0
        for i in range(len(pco)):
            for x in range(len(posScores)):
                if posScores[x] > pco[i]:
                    value +=1
            truePositivList.append(value)
            value = 0
        return truePositivList


    ## The TrueNegatives "TN" for a Probable Cut Off Value, ScoreList "negScores"
    def calculateTrueNegatives(self, pco, negScores):
        trueNegativesList = list()
        poc = pco
        negScores = negScores
        value = 0
        for i in range(len(pco)):
            for x in range(len(negScores)):
                if negScores[x] < poc[i]:
                    value +=1
            trueNegativesList.append(value)
            value = 0
        return trueNegativesList

    ## The TN-Rate "Specifity" i = length of the tnList, tnList = list of TN Values, totalNeg = total negative values
    def calculateTrueNegativRate(self, trueNegativList, totalNeg):
        trueNegativRateList = list()
        trueNegativList = trueNegativList
        totalNeg = totalNeg
        for i in range(len(trueNegativList)):
            value = float(trueNegativList[i]) / totalNeg
            trueNegativRateList.append(value)
        return trueNegativRateList

    ## The FP-Rate "1-TN Rate"
    def calculateFalsePositivRate(self, trueNegativRateList):
        falsePositivRateList = list()
        trueNegativRateList = trueNegativRateList
        for i in range(len(trueNegativRateList)):
            value = float(1 - trueNegativRateList[i])
            falsePositivRateList.append(value)
        return falsePositivRateList

    ## The TP-Rate "Sensivity" tpList = list of TP Values, totalPos = total positive values
    def truePositivRate(self, truePositivList, totalPos):
        truePositivRateList = list()
        truePositivList = truePositivList
        totalPos = totalPos
        for i in range(len(truePositivList)):
            value = float(truePositivList[i]) / totalPos
            truePositivRateList.append(value)
        return truePositivRateList

    ##########################################################################################################
    ##########################################################################################################
    ##                                                                                                      ##
    ##  To Calculate the FP-Rate List and the TP-Rate List:                                                 ##
    ##                                                                                                      ##
    ##  Use the probCutOff() to generate a List of Probable Cut Off Values                                  ##
    ##                                                                                                      ##
    ##  Iterate through the sorted Pos Scores using the truePositives() function and generate a List        ##
    ##  containing the values for each Probable Cut Off Value                                               ##
    ##                                                                                                      ##
    ##  Iterate through the sorted Neg Scores using the trueNegatives() function and generate a List        ##
    ##  containing the values for each Probable Cut Off Value                                               ##
    ##                                                                                                      ##
    ##  Iterate through the tnList using the tnRate() function and generate a List containing the values    ##
    ##  for each value                                                                                      ##
    ##                                                                                                      ##
    ##  Iterate through the tnRateList using the fpRate() function and generate a List containing           ##
    ##  the value for each value                                                                            ##
    ##                                                                                                      ##
    ##  Iterate through the tpList using the tpRate() function and generate a List containing the values    ##
    ##  for each value                                                                                      ##
    ##                                                                                                      ##
    ##########################################################################################################
    ##                                                                                                      ##
    ##                                                                                                      ##
    ##  Print the ROC Curve using the FP-Rate List on the X - Axis and the TP-Rate List on the Y - Axis     ##
    ##                                                                                                      ##
    ##                                                                                                      ##
    ##########################################################################################################
    ##########################################################################################################

    def setCutOff(self, probCutOffList):
        tmin = math.sqrt(self.falsePositivRateList[len(self.falsePositivRateList) - 1] ** 2 + (self.truePositivRateList[len(self.falsePositivRateList) - 1] - 1) ** 2)
        for i in range(len(probCutOffList)-1, -1, -1):
            t = math.sqrt(self.falsePositivRateList[i] ** 2 + (self.truePositivRateList[i] - 1) ** 2)
            if t < tmin:
                self.cutOffValue = probCutOffList[i]
                self.sensitivity = self.falsePositivRateList[i]
                self.specificity = self.truePositivRateList[i]



    def rocCurve(self, pcoStart, pcoEnd, pcoStep):
        print("generating the Data for the ROC Curve")
        pcoStart = pcoStart
        pcoEnd = pcoEnd
        pcoStep = pcoStep

        probCutOffList = self.probCutOff(pcoStart, pcoEnd, pcoStep)
        tpList = self.calculateTruePositives(probCutOffList, self.posScoresSorted)
        tnList = self.calculateTrueNegatives(probCutOffList, self.negScoresSorted)
        tnRateList = self.calculateTrueNegativRate(tnList, self.totalNegScores)
        self.falsePositivRateList = self.calculateFalsePositivRate(tnRateList)
        self.truePositivRateList = self.truePositivRate(tpList, self.totalPosScores)

        self.setCutOff(probCutOffList)

