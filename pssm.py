class Pssm:
    def __init__(self):
        self.frequencies =[]
        self.wmm = [{'A': 0.1, 'C': 0.1, 'G': 0.1, 'T': 0.1,},
                    {'A': 0.1, 'C': 0.1, 'G': 0.1, 'T': 0.1,},
                    {'A': 0.1, 'C': 0.1, 'G': 0.1, 'T': 0.1,},
                    {'A': 0.1, 'C': 0.1, 'G': 0.1, 'T': 0.1,}] # ight matrix model (rel freq)

        self.pssm = [{'A': -2.7, 'C': -1.5, 'G':-1.7000 , 'T': 1.7},
                     {'A': 1.8, 'C': -3.1, 'G': -4.9, 'T': -1.7},
                     {'A':0.1, 'C': -1.2, 'G': -1.1, 'T': 1.0},
                     {'A': 1.2, 'C': -1.0, 'G': -0.7, 'T': -1},
                     {'A': 1.0, 'C': -0.20, 'G': -1.1, 'T': -0.5},
                     {'A': -2.9, 'C': -2.2, 'G': -3.6, 'T': 1.8}] # position specific scoring matrix

    def calclateScore(self, motiv):
        score = 0
        for i in range(0, len(motiv)):
            score += self.pssm[i][motiv[i]]


#        for p in self.pssm:
#           score += p[motiv[0]]
#           motiv = motiv[1:]

        return score




if __name__  == "__main__":
    pssm = Pssm()
    print(Pssm().wmm[0]['A'])
    print pssm.calclateScore("TATAAT")