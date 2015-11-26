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

    ## Total Scores
    totalPosScores = len(pssm.posScores)
    totalNegScores = len(pssm.negScores)

    ## Sort the Scores
    posScoresSorted = sorted(pssm.posScores, reverse=True)
    negScoresSorted = sorted(pssm.negScores, reverse=True)

    ############################################### Functions ###############################################

    ## Generate a List of Probable Cut Off Values
    ## start = start Value, stop = stop Value, step = steps between the Values
    def probCutOff(start, stop, step):
        pronCutOffList = list()
        i = start
        while i < stop:
            pronCutOffList.append(i)
            i+= step
        return pronCutOffList

    cutOff = probCutOff(0.0, 8.0, 0.1)

    ## The TruePossitives "TP" for a Probable Cut Off Value,ScoreList "posScores" must be sorted
    def truePositives(pco, posScores):
        pco = pco
        value = 0
        for i in range(len(posScores)):
            if posScores[i] > pco:
                value +=1
        return value

    ## The TrueNegatives "TN" for a Probable Cut Off Value, ScoreList "negScores" must be sorted
    def trueNegatives(pco, negScores):
        poc = pco
        value = 0
        for i in range(len(negScores)):
            if negScores[i] > poc:
                value += 1
        value = len(negScores)-value
        return value

    ## The TN-Rate "Specifity" i = length of the tnList, tnList = list of TN Values, totalNeg = total negative values
    def tnRate(i, tnList, totalNeg):
        tnList = tnList
        value = (tnList[i] / totalNeg)*100
        return value

    ## The FP-Rate "1-TN Rate"
    def fpRate(i, tnRateList):
        tnRateList = tnRateList
        value = 100 - tnRateList[i]
        return value

    ## The TP-Rate "Sensivity" tpList = list of TP Values, totalPos = total positive values
    def tpRate(i, tpList, totalPos):
        tpList = tpList
        value = (tpList[i] / totalPos)*100
        return value

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
    ##  Iterate through the tnRateList using the fpPate() function and generate a List containing           ##
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