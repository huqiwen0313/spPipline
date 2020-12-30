import pandas as pd, numpy as np
from itertools import product


def create_jsoncodebook(infilepath, outfilepath, totalCycles=6, offCycles=2,
                        addEmptyBarcodes=True, uniColorAllowed=False,
                        firstCycleAnchor=False, anchorChannel=2):
    """
    Args:
        infilepath: a file containing gene names and barcodes (eg. 'SLC17A7_302301'). Each gene_barcode in a single-line
        outfilepath: path to the json codebook.
        totalCycles: number of total cycles to be decoded, including anchor and off-cycles.
        offCycles: number of off-cycles, i.e., the number of zeros in the barcode string.
        addEmptyBarcodes: If True, empty barcodes will be included (e.g. Empty_XXX where XXX is not in the input barcode list)
        uniColorAllowed: If False, enforces barcodes to have more than one color.
        Applicable only when addEmptyBarcodes==True.
        firstCycleAnchor: If True, the first cycle is assumed to be always on and
        the channel will be `anchorChannel`. Only effective if `addEmptyBarcodes is True.
    """

    # Generate Barcode List - each row a barcode
    barcodelist = np.array(list(product([0, 1, 2, 3], repeat=totalCycles)))

    # enforce the first round to be the anchor
    if firstCycleAnchor:
        barcodelist = barcodelist[barcodelist[:, 0] == anchorChannel]

    # enforce the number of off-cycles
    barcodelist = barcodelist[(barcodelist == 0).sum(axis=1) == offCycles,]

    # enforce the multi-color barcode
    barcodelist = pd.DataFrame(barcodelist)
    if not uniColorAllowed:
        if offCycles > 0:
            barcodelist = barcodelist.loc[barcodelist.nunique(axis=1) >= 3]
        else:
            barcodelist = barcodelist.loc[barcodelist.nunique(axis=1) >= 2]

    # making strings
    barcodelist = barcodelist.astype(str)
    barcode_str = barcodelist[0]
    for col in range(1, totalCycles):
        barcode_str = barcode_str + barcodelist[col]

    barcodelist = list(barcode_str)

    """ Adapted from Richard's codes"""
    barcode_dict = dict()
    # read Gene_Barcode file
    with open(infilepath) as file:
        for line in file:
            genename = line.strip('\n')
            barcode = genename.split('_')[1]
            barcode_dict[barcode] = genename
            if not barcode in barcodelist:
                raise ValueError("barcode {0} found in file {1} is not valid.".format(barcode, infilepath))

    # writing the json codebook
    codebook = open(outfilepath, 'w')
    codebook.write('{"version":"0.0.0","mappings":[')
    for bc in barcodelist:
        if bc in barcode_dict:
            genename = barcode_dict[bc]
        elif addEmptyBarcodes:
            genename = "Empty_{}".format(bc)

        codeword = convert_barcode_to_codeword(bc)
        codebook.write('{"codeword": ' + codeword + ', "target": "' + genename + '"}')
        if barcodelist.index(bc) < (len(barcodelist) - 1):
            codebook.write(',\n')

    codebook.write(']}')
    codebook.close()


def convert_barcode_to_codeword(barcode):
    """
    Convert string of 6 integers into starfish/json format codename
    barcode '1' = channel '0'
    barcode '2' = channel '1'
    barcode '3' = channel '2'

    302301 -> [{"c": 2, "r": 0, "v": 1.0}, {"c": 1, "r": 2, "v": 1.0}, {"c": 2, "r": 3, "v": 1.0}, {"c": 0, "r": 5, "v": 1.0}]
    """

    codeword = '['
    round = 0
    on = 0
    for i, character in enumerate(barcode):
        if character == '0':
            #             round = round + 1
            continue
        else:
            channel = str(int(character) - 1)
            codeword = codeword + '{"c": ' + channel + ', "r": ' + str(i) + ', "v": 1.0}, '

    codeword = codeword[0:-2] + ']'
    return codeword

