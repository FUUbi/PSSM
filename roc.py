from IPython.core.display import Math
import math


class Roc:
    def __init__(self, posScores, negScores):
        # Total Scores
        self.totalPosScores = len(posScores)
        self.totalNegScores = len(negScores)

        # Sort the Scores
        self.posScoresSorted = sorted(posScores, reverse=True)
        self.negScoresSorted = sorted(negScores, reverse=True)

        self.fpRateList = None
        self.tpRateList = None

        self.cutOffValue = None
        self.sensitivity = None
        self.specificity = None


    ## Generate a List of Probable Cut Off Values
    ## start = start Value, stop = stop Value, step = steps between the Values
    def __probCutOff(self, start, stop, step):
        pronCutOffList = list()
        i = start
        while i <= stop:
            i = round(i, 1)
            pronCutOffList.append(i)
            i+= step
        return pronCutOffList


    ## The TruePossitives "TP" for a Probable Cut Off Value,ScoreList "posScores" must be sorted
    def __truePositives(self, pco, posScores):
        tpList = list()
        pco = pco
        posScores = posScores
        value = 0
        for i in range(len(pco)):
            for x in range(len(posScores)):
                if posScores[x] > pco[i]:
                    value +=1
            tpList.append(value)
            value = 0
        return tpList


    ## The TrueNegatives "TN" for a Probable Cut Off Value, ScoreList "negScores" must be sorted
    def __trueNegatives(self, pco, negScores):
        tnList = list()
        poc = pco
        negScores = negScores
        value = 0
        for i in range(len(pco)):
            for x in range(len(negScores)):
                if negScores[x] < poc[i]:
                    value +=1
            tnList.append(value)
            value = 0
        return tnList

    ## The TN-Rate "Specifity" i = length of the tnList, tnList = list of TN Values, totalNeg = total negative values
    def __tnRate(self, tnList, totalNeg):
        tnRateList = list()
        tnList = tnList
        totalNeg = totalNeg
        for i in range(len(tnList)):
            value = float(tnList[i]) / totalNeg
            tnRateList.append(value)
        return tnRateList

    ## The FP-Rate "1-TN Rate"
    def __fpRate(self, tnRateList):
        fpRateList = list()
        tnRateList = tnRateList
        for i in range(len(tnRateList)):
            value = float(1 - tnRateList[i])
            fpRateList.append(value)
        return fpRateList

    ## The TP-Rate "Sensivity" tpList = list of TP Values, totalPos = total positive values
    def __tpRate(self, tpList, totalPos):
        tpRateList = list()
        tpList = tpList
        totalPos = totalPos
        for i in range(len(tpList)):
            value = float(tpList[i]) / totalPos
            tpRateList.append(value)
        return tpRateList

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

    def __setCutOff(self, probCutOffList):
        tmin = math.sqrt(self.fpRateList[len(self.fpRateList)-1]**2 + (self.tpRateList[len(self.fpRateList)-1] - 1)**2)
        for i in range(len(probCutOffList)-1, -1, -1):
            t = math.sqrt(self.fpRateList[i]**2 + (self.tpRateList[i] - 1)**2)
            if t< tmin:
                self.cutOffValue = probCutOffList[i]
                self.sensitivity = self.fpRateList[i]
                self.specificity = self.tpRateList[i]



    def rocCurve(self, pcoStart, pcoEnd, pcoStep):
        print("generatin the Data for the ROC Curve")
        print("If this takes more than 15 sec. please go to www.apple.com and by a daaaaaamn fucking MacBook Pro")
        print("or ... a sony vaio duo... daaaaam f**....and better buy it dont by it")
        pcoStart = pcoStart
        pcoEnd = pcoEnd
        pcoStep = pcoStep

        probCutOffList = self.__probCutOff(pcoStart, pcoEnd, pcoStep)
        tpList = self.__truePositives(probCutOffList, self.posScoresSorted)
        tnList = self.__trueNegatives(probCutOffList, self.negScoresSorted)
        tnRateList = self.__tnRate(tnList, self.totalNegScores)
        self.fpRateList = self.__fpRate(tnRateList)
        self.tpRateList = self.__tpRate(tpList, self.totalPosScores)

        self.__setCutOff(probCutOffList)

