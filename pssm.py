import Bio
from Bio import SeqIO
from Bio import motifs
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

        self.posScoresList = None
        self.negScoresList = None

    def calculatePosScore(self):
        print "calculating Pos Scores .... \t \t " ,
        self.posScoresList = list()
        for record in self.alignmentPos:
            pos = self.spliceSite.calculate(record.seq)
            neg = self.background.calculate(record.seq)
            self.posScoresList.append(pos - neg)
        print "done!"

    def calculateNegScore(self):
        print "calculating Neg Scores .... \t \t " ,
        self.negScoresList = list()
        for record in self.alignmentNeg:
            pos = self.spliceSite.calculate(record.seq)
            neg = self.background.calculate(record.seq)
            self.negScoresList.append(pos - neg)
        print "done!"

    def setSpliceSite(self, filePath):
        print "set Splice PSSM:  "
        self.spliceSite = self.calcutaltePssmSliceSite(filePath)

    def setBackground(self, filePath):
        print "set background PSSM: "
        self.background = self.calcutaltePssmBackground(filePath)

    def calcutaltePssmBackground(self, filePath):
        instance = list()
        for seq_record in SeqIO.parse(filePath, "fasta"):
            dna_seq = Seq(str(seq_record.seq).upper())
            instance.append(dna_seq)
        # wmm
        m = motifs.create(instance)
        m.weblogo("LogoBackground.png")
        pwm = m.counts.normalize(pseudocounts=0)
        pssm = pwm.log_odds()       # pssm
        print pssm
        return pssm

    def calcutaltePssmSliceSite(self, filePath):
        instance = list()
        for seq_record in SeqIO.parse(filePath, "fasta"):
            dna_seq = Seq(str(seq_record.seq).upper())
            instance.append(dna_seq)
        # wmm
        m = motifs.create(instance)
        m.weblogo("LogoSliceSite.png")
        pwm = m.counts.normalize(pseudocounts=0)
        pssm = pwm.log_odds()       # pssm
        print pssm
        return pssm

    def loadSeq(self, filePath):
        return Bio.AlignIO.read(filePath, "fasta", alphabet=IUPACUnambiguousDNA())

    def setAlignmentPos(self, filePath):
        print "loading AlignmentPos .... \t \t \t " ,
        self.alignmentPos = self.loadSeq(filePath)
        print "done!"

    def setAlignmentNeg(self, filePath):
        print "loading AlignmentNeg .... \t \t \t " ,
        self.alignmentNeg = self.loadSeq(filePath)
        print "done!"
