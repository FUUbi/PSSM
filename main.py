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

    cutOff = CutOff(pssm.posScoresList, pssm.negScoresList)
    cutOff.getCutOff(min(pssm.posScoresList), max(pssm.negScoresList))

    plt.hist(pssm.posScoresList, bins=100, color="blue", label="SpliceSites")
    plt.hist(pssm.negScoresList, bins=100, color="red", alpha=0.5, label="Background")
    plt.vlines(cutOff.value, 0, 4500, colors="green")
    plt.title("PSSM")

    plt.figure()
    plt.hist([s for s in pssm.posScoresList if s > cutOff.value], bins=100, color="blue", label="SpliceSites")
    plt.hist([s for s in pssm.negScoresList if s > cutOff.value], bins= 100, color="red", alpha=0.5, label="Background")
    plt.title("PSSM CUT OFF")

    print "Cut Off Value:" , cutOff.value
    print "Sensitivity: " + str(round(cutOff.sensitivityPercent, 2)) + "%"
    print "Specificity: " + str(round(cutOff.specificityPercent, 2)) + "%"

    print "Mcc: " + str(cutOff.mcc)

    plt.show()






