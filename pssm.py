from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq

# http://biopython-cn.readthedocs.org/en/latest/en/chr14.html

instance = []
for seq_record in SeqIO.parse("data/test_neg.txt", "fasta"):
    dna_seq = Seq(str(seq_record.seq).upper())
    instance.append(dna_seq)

# wmm
m = motifs.create(instance)
pwm = m.pwm
print(pwm)

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

    def calclateScore(self, motiv):
        score = 0
        for i in range(0, len(motiv)):
            score += self.pssm[i][motiv[i]]


#        for p in self.pssm:
#           score += p[motiv[0]]
#           motiv = motiv[1:]

        return score

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



if __name__  == "__main__":
    pssm = Pssm()
    pssm.setSpliceSite("data/train_pos.txt")
    pssm.setBackground("data/train_neg.txt")


    print(pssm.background)
    print(pssm.spliceSite)

