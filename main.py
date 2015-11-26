from pssm import *
from roc import *

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

    roc = Roc(pssm.posScores, pssm.negScores)
    roc.rocCurve(-1.0, 4.0, 0.1)

    print(roc.cutOffValue)
    print(roc.sensitivity)
    print(roc.specificity)
    plt.plot(roc.fpRateList, roc.tpRateList, color="purple")
    plt.show()





