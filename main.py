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

    roc2 = Roc2(pssm.posScoresList, pssm.negScoresList)
    roc2.setCutOff(-1, 4)
    #cutOff = CutOff(pssm.posScoresList, pssm.negScoresList)
    #cutOffValue = cutOff.getCutOff(-1, 4)

    plt.hist(pssm.posScoresList, bins=100, color="blue", label="SpliceSites")
    plt.hist(pssm.negScoresList, bins=100, color="red", alpha=0.5, label="Background")
    plt.vlines(roc2.cuOff.value,0,4500, colors="green")
    #plt.vlines(cutOffValue, 0, 4500, colors="green")
    # plt.vlines(roc.cutOffValue,0,4500, colors="green")
    plt.title("PSSM")

    plt.figure()
    plt.hist([s for s in pssm.negScoresList if s > cutoffValue], bins= 100, color="red")
    plt.hist([s for s in pssm.posScoresList if s > cutoffValue], bins=100, color="blue")
    #plt.hist([s for s in pssm.posScoresList if s > cutOffValue], bins=100, color="blue", label="SpliceSites")
    #plt.hist([s for s in pssm.negScoresList if s > cutOffValue], bins= 100, color="red", alpha=0.5, label="Background")
    plt.title("PSSM CUT OFF")

    # plt.figure()
    # plt.plot(roc.falsePositiveRateList, roc.truePositiveRateList, color="purple")
    # plt.title("ROC")

    print "Cut Off Value:" , roc2.cuOff.value
    # print "Cut Off Value:" , cutOffValue
    # print "Sensitivity:" , roc.sensitivity
    #print "Sensitivity: " + str(round(cutOff.sensitivityPercent, 2)) + "%"
    # print "Specificity:" , roc.specificity
    #print "Specificity: " + str(round(cutOff.specificityPercent, 2)) + "%"

    plt.show()






