import math

import scipy.optimize as so

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

        #http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html#scipy.optimize.differential_evolution
        self.mcc = so.differential_evolution(self.calcMcc,  [(probableCutOffStart, probableCutOffEnd)]).x[0]

        #### more random optimation, also works.....
        # func = lambda x: self.calcMcc(x)
        # x0 =[1.] #Initial guess.
        # minimizer_kwargs = {"method": "BFGS"}
        # ret = basinhopping(func, x0, minimizer_kwargs=minimizer_kwargs, niter=100)
        # self.mcc = ret.x



        self.calcSensitivity()
        self.calcSpecificity()
        print "done!"


    def calcMcc(self, x):
        truePositive = self.calculatePositive(self.posScoresList, x)
        falsePositive = self.calculatePositive(self.negScoresList, x)

        trueNegative = self.totalNegScores - falsePositive
        falseNegative = self.totalPosScores - truePositive

        return -float(truePositive * trueNegative - falsePositive * falseNegative) / \
                 (math.sqrt((truePositive + falsePositive) * (trueNegative + falsePositive) *
                            (truePositive + falseNegative) * (trueNegative + falseNegative))
                  )
