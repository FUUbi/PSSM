from pssm import *
from cutoff import *

if __name__  == "__main__":
    pssm = Pssm()
    pssm.setSpliceSite("data/train_pos.txt")
    pssm.setBackground("data/train_neg.txt")
    pssm.setAlignmentPos("data/test_pos.txt")
    pssm.setAlignmentNeg("data/test_neg.txt")

    pssm.calculateNegScore()
    pssm.calculatePosScore()

#    roc = Roc(pssm.posScoresList, pssm.negScoresList)
#    roc.rocCurve(-1.0, 4.0, 0.1)

    cutoff = CutOff(pssm.posScoresList, pssm.negScoresList)
    cutoffValue = cutoff.getCutOff(-1, 4)

    plt.hist(pssm.posScoresList, bins=100, color="blue", label="SpliceSites")
    plt.hist(pssm.negScoresList, bins=100, color="red", alpha=0.5, label="Background")
    plt.vlines(cutoffValue,0,4500, colors="green")
    # plt.vlines(roc.cutOffValue,0,4500, colors="green")
    plt.title("PSSM")

    plt.figure()
    plt.hist([s for s in pssm.negScoresList if s > .5], bins= 100, color="red")
    plt.hist([s for s in pssm.posScoresList if s > .5], bins=100, color="blue")
    plt.title("PSSM CUT OFF")

    # plt.figure()
    # plt.plot(roc.falsePositiveRateList, roc.truePositiveRateList, color="purple")
    # plt.title("ROC")

    print "Cut Off Value:" , cutoffValue
    # print "Sensitivity:" , roc.sensitivity
    # print "Specificity:" , roc.specificity

    plt.show()






