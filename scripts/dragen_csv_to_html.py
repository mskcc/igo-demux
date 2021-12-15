import pandas
import sys

# Combines the DRAGEN Demultiplex_Stats.csv and Top_Unknown_Barcodes.csv into one .html file
if __name__ == "__main__":
    #Usage: python dragen_csv_to_html.py [dragen_demux_dir] [output_file_name]
    demultiplex_stats = sys.argv[1] + "Demultiplex_Stats.csv"
    top_unknown_barcodes = sys.argv[1] + "Top_Unknown_Barcodes.csv" 
    write_to_file = sys.argv[2]

    print("Converting DRAGEN {} .csv files to .html file".format(demultiplex_stats))

    demux_stats_csv = pandas.read_csv(demultiplex_stats)
    top_unknown_barcodes_csv = pandas.read_csv(top_unknown_barcodes)

    # not identical to original laneBarcode.html but still mostly readable
    # TODO formate two tables in the html with different column headers
    final_df = pandas.concat([demux_stats_csv, top_unknown_barcodes_csv])

    final_df.to_html(write_to_file)