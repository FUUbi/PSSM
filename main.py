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

    roc = Roc(pssm.posScores, pssm.negScores)
    roc.rocCurve(-1.0, 4.0, 0.1)

    plt.hist(pssm.negScores, bins= 100, color="red")
    plt.hist(pssm.posScores, bins=100, color="blue")
    plt.vlines(roc.cutOffValue,0,4500, colors="green")
    plt.title("PSSM")

    plt.figure()
    plt.hist([s for s in pssm.negScores if s > .5], bins= 100, color="red")
    plt.hist([s for s in pssm.posScores if s > .5], bins=100, color="blue")
    plt.title("PSSM CUT OFF")

    plt.figure()
    plt.plot(roc.falsePositivRateList, roc.truePositivRateList, color="purple")
    plt.title("ROC")

    print(roc.cutOffValue)
    print(roc.sensitivity)
    print(roc.specificity)

    plt.show()






