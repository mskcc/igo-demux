import pandas
import sys

from pandas.io import html

if __name__ == "__main__":
    demultiplex_stats = sys.argv[1] + "Demultiplex_Stats.csv"
    top_unknown_barcodes = sys.argv[1] + "Top_Unknown_Barcodes.csv" 
    write_to_file = sys.argv[2]

    print("Converting DRAGEN {} .csv files to .html file".format(demultiplex_stats))

    demux_stats_csv = pandas.read_csv(demultiplex_stats)
    top_unknown_barcodes_csv = pandas.read_csv(top_unknown_barcodes)

    # not identical to original laneBarcode.html but still mostly readable
    final_df = pandas.concat([demux_stats_csv, top_unknown_barcodes_csv])

    final_df.to_html(write_to_file)