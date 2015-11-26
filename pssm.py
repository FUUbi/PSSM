from Bio import SeqIO
from Bio import motifs
import Bio
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
