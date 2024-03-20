def folder_name2title(folder_name):
    cutoff=float(folder_name.split("_")[1])
    return "$\Delta_{\mathrm{cut}}="+str(cutoff)+"$"

def trotter_name2title(trotter_name):
    divider, divided = trotter_name[7:].split("-")
    if divider != "":
        return "Trotter {}/{}".format(divided,divider)
    else:
        return "Trotter {}".format(divided)