import pandas
import sys
from collections import OrderedDict

# Combines the DRAGEN Demultiplex_Stats.csv and Top_Unknown_Barcodes.csv into one .html file
def dragen_csv_to_html(dragen_demux_dir, output_file_name):
    demultiplex_stats = dragen_demux_dir + "Demultiplex_Stats.csv"
    top_unknown_barcodes = dragen_demux_dir + "Top_Unknown_Barcodes.csv" 
    write_to_file = output_file_name

    print("Converting DRAGEN {} .csv files to .html file".format(demultiplex_stats))

    demux_stats_csv = pandas.read_csv(demultiplex_stats)
    # convert int dtype to float in order to add commas to reads number
    demux_stats_csv_convert = demux_stats_csv.astype({"# Reads": 'float64', "# Perfect Index Reads": 'float64', "# One Mismatch Index Reads": 'float64'})

    top_unknown_barcodes_csv = pandas.read_csv(top_unknown_barcodes)
    # convert int dtype to float in order to add commas to reads number
    top_unknown_barcodes_csv_covert = top_unknown_barcodes_csv.astype({"# Reads": 'float64'})
    # seperate top unknown barcode according to lane number and store in orderedDict named df_by_lanes
    lane_number = max(top_unknown_barcodes_csv_covert["Lane"])
    df_by_lanes = OrderedDict()
    for i in range(1, lane_number + 1):
        df_name = "top_unknown_lane" + str(i)
        df_by_lanes[df_name] = top_unknown_barcodes_csv_covert.loc[top_unknown_barcodes_csv_covert["Lane"] == i]
    
    # format two tables in the html with different column headers
    with open(write_to_file, 'w') as _file:
        _file.write("<h2>Lane Summary<h2>" + demux_stats_csv_convert.to_html(index = False, float_format =  '{:,.0f}'.format) + "\n<h2>Top Unknown Barcodes<h2>\n" + "<table>\n" )
        for value in df_by_lanes.values():
            _file.write("<td>" + value.to_html(index = False, float_format =  '{:,.0f}'.format) + "</td>")
        _file.write("\n</table>")


if __name__ == "__main__":
    #Usage: python dragen_csv_to_html.py [dragen_demux_dir] [output_file_name]
    dragen_csv_to_html(sys.argv[1], sys.argv[2])