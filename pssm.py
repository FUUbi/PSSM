from Bio import SeqIO
from Bio import motifs
import Bio
from Bio.Alphabet import IUPAC
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA
from Bio.Seq import Seq
import matplotlib.pyplot as plt

# http://biopython-cn.readthedocs.org/en/latest/en/chr14.html
class Pssm:
    def __init__(self):
        self.spliceSite = None
        self.background = None

        self.alignmentPos = None
        self.alignmentNeg = None

        self.posScores = None
        self.negScores = None

    def calculatePosScore(self):
        self.posScores = []
        for record in self.alignmentPos:
            self.posScores.append(self.spliceSite.calculate(record.seq))
        print "calculated Pos Scores"

    def calculateNegScore(self):
        self.negScores = []
        for record in self.alignmentNeg:
            self.negScores.append(self.background.calculate(record.seq))
        print "calculated Neg Scores"

    def setSpliceSite(self, filePath):
        self.spliceSite = self.calcutaltePssm(filePath)
        print "set Splice PSSM  "

    def setBackground(self, filePath):
        self.background = self.calcutaltePssm(filePath)
        print "set background PSSM"

    def calcutaltePssm(self, filePath):
        instance = []
        for seq_record in SeqIO.parse(filePath, "fasta"):
            dna_seq = Seq(str(seq_record.seq).upper())
            instance.append(dna_seq)
        # wmm
        m = motifs.create(instance)

        pwm = m.counts.normalize(pseudocounts=0.5)
        return pwm.log_odds()       # pssm

    def loadSeq(self, filePath):
        return Bio.AlignIO.read(filePath, "fasta", alphabet=IUPACUnambiguousDNA())

    def setAlignmentPos(self, filePath):
        self.alignmentPos = self.loadSeq(filePath)

        print "loaded AlignmentPos"

    def setAlignmentNeg(self, filePath):
        self.alignmentNeg = self.loadSeq(filePath)
        print "loaded AlignmentNeg"

if __name__  == "__main__":
    pssm = Pssm()
    pssm.setSpliceSite("data/train_pos.txt")
    pssm.setBackground("data/train_neg.txt")
    pssm.setAlignmentPos("data/test_pos.txt")
    pssm.setAlignmentNeg("data/test_neg.txt")

    pssm.calculateNegScore()
    pssm.calculatePosScore()

    plt.hist(pssm.negScores, bins= 100, color="red")
    plt.hist(pssm.posScores, bins=100, color="blue")

    plt.show()



    #########################################################################################################
    #########################################################################################################
    #########################################################################################################
    # negScoresSorted = list()
    # posScoresSorted = list()
    #
    # negScoresSorted = (3.19113753652, 1.76896238158, 1.10029167012, 0.25975834748, -1.44403621830, -2.18898610220, -2.60898178181, -3.13453393946, -3.22828617301, -4.77617608561)
    #
    # posScoresSorted = (6.84946860220, 5.21223598584, 5.05555992006, 4.69326858755, 4.16360951066, 3.68749250459, 3.51298798280, 2.45355113003, 1.46175619017, 0.363484753632)
    #
    # print(negScoresSorted)
    # print(posScoresSorted)
    #
    # totalNegScores = len(negScoresSorted)
    # totalPosScores = len(posScoresSorted)

    # Total Scores
    totalPosScores = len(pssm.posScores)
    totalNegScores = len(pssm.negScores)

    # Sort the Scores
    posScoresSorted = sorted(pssm.posScores, reverse=True)
    negScoresSorted = sorted(pssm.negScores, reverse=True)

    ############################################### Functions ###############################################

    ## Generate a List of Probable Cut Off Values
    ## start = start Value, stop = stop Value, step = steps between the Values
    def probCutOff(start, stop, step):
        pronCutOffList = list()
        i = start
        while i <= stop:
            i = round(i, 1)
            pronCutOffList.append(i)
            i+= step
        return pronCutOffList


    ## The TruePossitives "TP" for a Probable Cut Off Value,ScoreList "posScores" must be sorted
    def truePositives(pco, posScores):
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
    def trueNegatives(pco, negScores):
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
    def tnRate(tnList, totalNeg):
        tnRateList = list()
        tnList = tnList
        totalNeg = totalNeg
        for i in range(len(tnList)):
            value = float(tnList[i]) / totalNeg
            tnRateList.append(value)
        return tnRateList

    ## The FP-Rate "1-TN Rate"
    def fpRate(tnRateList):
        fpRateList = list()
        tnRateList = tnRateList
        for i in range(len(tnRateList)):
            value = float(1 - tnRateList[i])
            fpRateList.append(value)
        return fpRateList

    ## The TP-Rate "Sensivity" tpList = list of TP Values, totalPos = total positive values
    def tpRate(tpList, totalPos):
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


    def rocCurve(pcoStart, pcoEnd, pcoStep, sortedPos, sortedNeg):
        print("generatin the Data for the ROC Curve")
        print("If this takes more than 15 sec. please go to www.apple.com and by a daaaaaamn fucking MacBook Pro")

        pcoStart = pcoStart
        pcoEnd = pcoEnd
        pcoStep = pcoStep
        posScoresSorted = sortedPos
        negScoresSorted = sortedNeg

        probCutOffList = probCutOff(pcoStart, pcoEnd, pcoStep)
        tpList = truePositives(probCutOffList, posScoresSorted)
        tnList = trueNegatives(probCutOffList, negScoresSorted)
        tnRateList = tnRate(tnList, totalNegScores)
        fpRateList = fpRate(tnRateList)
        tpRateList = tpRate(tpList, totalPosScores)

        # print(probCutOffList)
        #
        # print(tpList)
        #
        # print(tnList)
        #
        # print(tnRateList)
        #
        # print(fpRateList)
        #
        # print(tpRateList)

        fp = fpRateList
        tp = tpRateList

        plt.plot(fp, tp, color="purple")
        plt.show()






    rocCurve(-1.0, 4.0, 0.1, posScoresSorted, negScoresSorted)


