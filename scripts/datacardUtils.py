from decimal import Decimal
def GetNDecimalDigits(num):
    d = Decimal(str(num))
    positiveExponent = abs(d.as_tuple().exponent)
    return positiveExponent

def RoundToN(x, n, returnFloat=False):
    # if n < 1:
    #    raise ValueError("can't round to less than 1 sig digit!")
    # # number of digits given by n
    # return "%.*e" % (n-1, x)
    if isinstance(x, float):
        rounded = round(x, n)
        # return "{0:.{1}f}".format(rounded, n)
        if returnFloat:
            return rounded
        # return str(rounded)
        return "{0:.{1}f}".format(rounded, n)
    else:
        return str(x)

def FormatToNDigits(num, n):
    if isinstance(num, str) or isinstance(num, int):
        return str(num)
    rounded = RoundToN(num, n, returnFloat=True)
    return "{0:.{1}f}".format(rounded, n)

def RoundToNSigFigs(num, n=1):
    if isinstance(num, str):
        return num, -1
    nDecDigits = GetNDecimalDigits(num)
    d = Decimal(str(num))
    nDigits = len(d.normalize().as_tuple().digits)
    digitsBeforeDecimal = nDigits - nDecDigits
    d = round(d, n-digitsBeforeDecimal)
    # dStr = str(np.format_float_positional(float(d), trim='-'))
    dStr = str(int(d))
    if d < 10:
        dStr = "{0:.{1}f}".format(float(d), n-digitsBeforeDecimal)
    # print("DEBUG: RoundToNSigFigs for num={} --> {}; n-digitsBeforeDecimal = {}-{} = {}".format(num, dStr, n, digitsBeforeDecimal, n-digitsBeforeDecimal))
    # but now, digits before decimal could have changed, e.g., 0.99 --> 1.0
    if "." in dStr:
        d = Decimal(dStr)
        nDecDigitsUpdated = GetNDecimalDigits(float(d))
        if dStr[-1] == "0":
            nDecDigitsUpdated += 1 # handle case with trailing zero
        dUpdated = Decimal(dStr)
        # nDigitsUpdated = len(dUpdated.as_tuple().digits) if "." in dStr else len(d.normalize().as_tuple().digits)
        nDigitsUpdated = len(dUpdated.normalize().as_tuple().digits)
        # print("DEBUG: RoundToNSigFigs for num={} updated: d={}, nDecDigitsUpdated={}, dUpdated={}, nDigitsUpdated={}".format(num, d, nDecDigitsUpdated, dUpdated, nDigitsUpdated))
        digitsBeforeDecimalUpdated = nDigitsUpdated - nDecDigitsUpdated
        if digitsBeforeDecimalUpdated > digitsBeforeDecimal:
            dStr = str(int(d))
        # print("DEBUG: RoundToNSigFigs for num={} --> {}; [updated] n-digitsBeforeDecimalUpdated = {}-{} = {}".format(num, dStr, n, digitsBeforeDecimalUpdated, n-digitsBeforeDecimalUpdated))
    return dStr, n-digitsBeforeDecimal

def GetTableEntryStr(evts, errStatUp="-", errStatDown="-", errSyst=0, addDecimalsUntilNonzero=False, latex=False):
    if evts == "-":
        return evts
    #print("DEBUG: [1] GetTableEntryStr() RoundToNSigFigs for errStatUp={}".format(errStatUp))
    errStatUpR, digitsAwayFromDecimal = RoundToNSigFigs(errStatUp)
    #print("DEBUG: GetTableEntryStr() AFTER RoundToNSigFigs for errStatUp={}: got errStatUpR={}, digitsAwayFromDecimal={}".format(errStatUp, errStatUpR, digitsAwayFromDecimal))
    errStatDownR = str(float(round(Decimal(errStatDown), digitsAwayFromDecimal))) if not isinstance(errStatDown, str) else errStatDown
    # print("DEBUG: GetTableEntryStr() Now for errStatDown={}: got errStatDownR={}, digitsAwayFromDecimal={}".format(errStatDown, errStatDownR, digitsAwayFromDecimal))
    evtsR = str(float(round(Decimal(evts), digitsAwayFromDecimal)))
    if float(evtsR) > 1 and evtsR.endswith(".0") and len(evtsR) > len(errStatUpR):
        evtsR = evtsR[:-2]
    elif errStatUpR != "-" and len(evtsR) != len(errStatUpR):
        #print("DEBUG: GetTableEntryStr() reformat evtsR to {} digits: evtsR={} --> {}".format(GetNDecimalDigits(errStatUpR), evtsR, "{0:.{1}f}".format(float(evtsR), GetNDecimalDigits(errStatUpR))))
        evtsR = "{0:.{1}f}".format(float(evtsR), GetNDecimalDigits(errStatUpR))
    #print("DEBUG: GetTableEntryStr() Now for evts={}: got evtsR={}, digitsAwayFromDecimal={}".format(evts, evtsR, digitsAwayFromDecimal))
    # # rounding
    # errStatUpR = RoundToN(errStatUp, 2)
    # if GetNDecimalDigits(errStatUpR) > 1:
    #     errStatDownR = FormatToNDigits(errStatDown, 2)
    #     evtsR = FormatToNDigits(evts, 2)
    #     errStatUpR = FormatToNDigits(errStatUp, 2)
    # else:
    #     errStatDownR = FormatToNDigits(errStatDown, 1)
    #     evtsR = FormatToNDigits(evts, 1)
    #     errStatUpR = FormatToNDigits(errStatUp, 1)

    # suppress
    suppressed = False
    if errStatUpR != "-" and float(errStatUpR) < 0.01:
        errStatUpR = "0.0"
        suppressed = True
    if errStatDownR != "-" and float(errStatDownR) < 0.01:
        errStatDownR = "0.0"
        suppressed = True
    if suppressed:
        evtsR, _ = RoundToNSigFigs(evts, n=2)
        if float(evtsR) < 0.01:
            evtsR, _ = RoundToNSigFigs(evts, n=1)
    # add additional decimal place if it's zero after rounding
    if addDecimalsUntilNonzero and float(evtsR) == 0.0:
        # print("DEBUG: adding decimals until nonzero for evts={}; evtsR={}; errStatUp={}, errStatDown={}".format(evts, evtsR, errStatUp, errStatDown))
        nDigits = 2
        while float(evtsR) == 0.0 and nDigits < 5:
            nDigits += 1
            evtsR = RoundToN(evts, nDigits)
            # print("\tDEBUG: rounded to n={}, evtsR={}".format(nDigits, evtsR))
        # errStatUpR = RoundToN(errStatUp, nDigits)
        # errStatDownR = RoundToN(errStatDown, nDigits)
        # print("\tDEBUG: for evts={}; ended up with evtsR={}; errStatUpR={}, errStatDownR={}".format(evts, evtsR, errStatUpR, errStatDownR))
    # handle cases where we don't specify stat or syst
    #print(errStatUp,errSyst)
    if errStatUp == "-":
        return evtsR
    elif errSyst == 0:
        #print("errsyst = 0")
        if errStatUp == errStatDown:
            if not latex:
                return evtsR + " +/- " + errStatUpR if float(errStatUpR) != 0 else evtsR
            else:
                return evtsR + " \\pm " + errStatUpR if float(errStatUpR) != 0 else evtsR
        else:
            if not latex:
                return evtsR + " + " + errStatUpR + " - " + errStatDownR
            else:
                return (
                    evtsR
                    + "^{+"
                    + errStatUpR
                    + "}_{-"
                    + errStatDownR
                    + "}"
                )
    else:
        #print("errsyst != 0")
        # errSystR = FormatToNDigits(errSyst, 2)
        errSystR = str(float(round(Decimal(errSyst), digitsAwayFromDecimal))) if not isinstance(errSyst, str) else errSyst
        # suppress or remove trailing zero
        if errSystR != "-" and float(errSystR) < 0.01:
            errSystR = "0.0"
        elif float(errSystR) > 1 and errSystR.endswith(".0") and len(errSystR) > len(evtsR):
            errSystR = errSystR[:-2]
        elif errStatUpR != "-" and len(errSystR) != len(errStatUpR):
            # print("DEBUG: reformat errSystR to {} digits: errSystR={} --> {}".format(GetNDecimalDigits(errStatUpR), errSystR, "{0:.{1}f}".format(float(errSystR), GetNDecimalDigits(errStatUpR))))
            errSystR = "{0:.{1}f}".format(float(errSystR), GetNDecimalDigits(errStatUpR))
        if errStatUp == errStatDown:
            #print("errstatup = errstatdown")
            if not latex:
                return evtsR + " +/- " + errStatUpR + " +/- " + errSystR
            else:
                #print("will return "+ evtsR + " \\pm " + errStatUpR + " \\pm " + errSystR)
                return (
                    evtsR + " \\pm " + errStatUpR + " \\pm " + errSystR
                )
        else:
            return (
                evtsR
                + "^{+"
                + errStatUpR
                + "}_{-"
                + errStatDownR
                + "} \\pm "
                + errSystR
            )