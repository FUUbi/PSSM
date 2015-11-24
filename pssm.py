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
pwm = m.counts.normalize(pseudocounts=0.5)

# You can also directly access columns of the counts matrix
m.counts[:,3]

# You can access these counts as a dictionary:
m.counts['A']

# pssm
pssm =pwm.log_odds()
