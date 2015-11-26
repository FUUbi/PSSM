from Bio import SeqIO
from Bio import motifs
import Bio
from Bio.Alphabet import IUPAC
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA
from Bio.Seq import Seq
import matplotlib.pyplot as plt

# http://biopython-cn.readthedocs.org/en/latest/en/chr14.html


instance = []
for seq_record in SeqIO.parse("data/test_neg.txt", "fasta"):
    dna_seq = Seq(str(seq_record.seq).upper())
    instance.append(dna_seq)


# wmm
m = motifs.create(instance)
pwm = m.pwm
print(pwm)
print(dir(m))

# You can also directly access columns of the counts matrix
print  m.counts[:,3]

# You can access these counts as a dictionary:
print  m.counts['A']

# pssm
pssm = pwm.log_odds()
print pssm['A',3]

print pssm
class Pssm:
    def __init__(self):
        self.spliceSite = None
        self.background = None

        self.seqPos = None
        self.seqNeg = None

        self.posScores = None
        self.negScores = None

    def calculatePosScore(self):
        self.posScores = []

        for seq in self.seqPos:
#          self.posScores.append(self.spliceSite.calculate(seq))
            score = 0
            for i in range(0, len(seq)):
                score += self.spliceSite[seq[i], i]
                self.posScores.append(score)

    def calculateNegScore(self):
        self.negScores = []
        for seq in self.seqNeg:
            #self.negScores.append(self.background.calculate(seq))
            score = 0
            for i in range(0, len(seq)):
                score += self.spliceSite[seq[i], i]
            self.negScores.append(score)

    def setSpliceSite(self, filePath):
        self.spliceSite = self.calcutaltePssm(filePath)

    def setBackground(self, filePath):
        self.background = self.calcutaltePssm(filePath)

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
        instance = []
        for seq_record in SeqIO.parse(filePath, "fasta"):
            dna_seq = Seq(str(seq_record.seq).upper(), IUPACUnambiguousDNA())
            instance.append(dna_seq)

        return instance

    def setSeqPos(self, filePath):
        self.seqPos = self.loadSeq(filePath)

    def setSeqNeg(self, filePath):
        self.seqNeg = self.loadSeq(filePath)

if __name__  == "__main__":
    pssm = Pssm()
    pssm.setSpliceSite("data/train_pos.txt")
    pssm.setBackground("data/train_neg.txt")
    pssm.setSeqPos("data/test_pos.txt")
    pssm.setSeqNeg("data/test_neg.txt")

    print(pssm.background)
    print(pssm.spliceSite)

    pssm.calculateNegScore()
    pssm.calculatePosScore()
    for i in range(0, 10):
        print("neg" + str(pssm.negScores[i]) + "   pos " + str(pssm.posScores[i]))
    plt.hist(pssm.negScores, bins= 100, color="red")
    plt.hist(pssm.posScores, bins= 100, color="blue")

    #plt.show()

    alignment = Bio.AlignIO.read("data/test_pos.txt", "fasta")
    for record in alignment :
        print record.seq, record.id
        print dir(record)