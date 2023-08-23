# Code that re-writes Illumina HTML Reports wich have barcodes so the barcodes have tooltips with their names 
# i.e. an HTML report with DNA barcode "GCGCAAGC-TCACGCCG" has tooltip "UDI0042"
# if the sequncer has reversed i5, add reversed as third input
# Usage: python3 barcodelookup.py laneBarcode.html Barcodes.json reversed

import re
import json
import sys


def reverse_barcode(barcode):
    # given a barcode such as GGGGGGGG-AGATCTCG, reverse complementary the part after - and return
    barcode_lst = barcode.split('-')
    if len(barcode_lst) != 2:
        return barcode
    else:
        barcode_reversed = []
        for letter in barcode_lst[1]:
            if letter == "A":
                barcode_reversed.insert(0, "T")
            elif letter == "T":
                barcode_reversed.insert(0, "A")
            elif letter == "C":
                barcode_reversed.insert(0, "G")
            elif letter == "G":
                barcode_reversed.insert(0, "C")
        barcode_reversed = "".join(barcode_reversed)
        barcode_reversed = "-".join([barcode_lst[0], barcode_reversed])
    return barcode_reversed

if (len(sys.argv) > 1):
    htmlfile = sys.argv[1]

if __name__ == "__main__":
    jsonfile = sys.argv[2] # The dictionary of known barcodes
    newhtmlstring = htmlfile.replace(".html", "_.html") # The name of the new HTML file to write to
    newhtml = open(newhtmlstring,'w')
        
    # build barcode dictionary from Barcodes.json file
    with open(jsonfile) as f:
        jsondata = json.load(f)
    bcdict = {}
    if sys.argv[3] == "reversed":
        for key in jsondata:
            key['INDEXTAG'] = reverse_barcode(key['INDEXTAG'])
            if key['INDEXTAG'] in bcdict:
                bcdict[key['INDEXTAG']].append(key['INDEXID'])
            else:
                bcdict[key['INDEXTAG']] = [key['INDEXID']]

    else:
        for key in jsondata:
            if key['INDEXTAG'] in bcdict:
                bcdict[key['INDEXTAG']].append(key['INDEXID'])
            else:
                bcdict[key['INDEXTAG']] = [key['INDEXID']]
    #taken from https://www.w3schools.com/css/css_tooltip.asp
    styletext = """
    <style>
    /* Tooltip container */
    .tooltip {
        position: relative;
        display: inline-block;
        border-bottom: 1px dotted black; /* If you want dots under the hoverable text */
    }
    /* Tooltip text */
    .tooltip .tooltiptext {
        visibility: hidden;
        width: 120px;
        background-color: black;
        color: #fff;
        text-align: center;
        padding: 5px 0;
        border-radius: 6px;
    
        /* Position the tooltip text - see examples below! */
        position: absolute;
        z-index: 1;
    }
    /* Show the tooltip text when you mouse over the tooltip container */
    .tooltip:hover .tooltiptext {
        visibility: visible;
    }
    </style>
    """

    f = open(htmlfile).read().splitlines()
    # The regex to look for all DNA barcodes in the Illumina HTML reports
    reglist = [re.findall(r'\<td\>([ATGC]+\-?[ATGC]+)\<\/td\>', line) for line in f]
    assert len(reglist) == len(f)
    for idx,line in enumerate(reglist):
        if idx == 2:
            print(f[idx], file=newhtml)
            print(styletext, file=newhtml) # add the css for the tooltip
        elif 'hide barcodes' in f[idx]: # remove 'hide barcodes' link because not all HTML reports are copied, link is broken
            pass 
        else:
            if len(line) == 1:  # one barcode expected per line in input HTML
                bckey = line[0].replace('+','-')
                try:
                    tooltipkey = ",".join(bcdict[bckey])
                    newstring = '<td><div class="tooltip">'+line[0]+'<span class="tooltiptext">'+tooltipkey+'</span></div></td>'
                    print(newstring, file=newhtml)
                except KeyError:
                    print(f[idx], file=newhtml)
            else:
                print(f[idx], file=newhtml)
    newhtml.close()