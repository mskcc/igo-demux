from cmath import nan
import pandas
import sys
from collections import OrderedDict

""" 
Combines the DRAGEN Reports/Demultiplex_Stats.csv and Reports/Top_Unknown_Barcodes.csv into one laneBarcode.html file imitating the output from bcl2fastq
"""
def build_lane_summary_html(demux_reports_dir, write_to_file):
    demultiplex_stats = demux_reports_dir + "/Demultiplex_Stats.csv"
    top_unknown_barcodes = demux_reports_dir + "/Top_Unknown_Barcodes.csv" 

    print("Converting DRAGEN {} and {} .csv files to .html file".format(demultiplex_stats, top_unknown_barcodes))

    demux_stats_csv = pandas.read_csv(demultiplex_stats)
    # convert int dtype to float in order to add commas to reads number
    demux_stats_csv_convert = demux_stats_csv.astype({"# Reads": 'float64', "# Perfect Index Reads": 'float64', "# One Mismatch Index Reads": 'float64'})
    demux_stats_csv_convert["% Reads"] = demux_stats_csv_convert["% Reads"] * 100
    top_unknown_barcodes_csv = pandas.read_csv(top_unknown_barcodes)
    # convert int dtype to float in order to add commas to reads number
    top_unknown_barcodes_csv_covert = top_unknown_barcodes_csv.astype({"# Reads": 'float64'})
    # seperate top unknown barcode according to lane number and store in orderedDict named df_by_lanes
    lane_number = max(top_unknown_barcodes_csv_covert["Lane"])
    df_by_lanes = OrderedDict()
    for i in range(1, lane_number + 1):
        df_name = "top_unknown_lane" + str(i)
        df_by_lanes[df_name] = top_unknown_barcodes_csv_covert.loc[top_unknown_barcodes_csv_covert["Lane"] == i]
        if not df_by_lanes[df_name]["index2"].isnull().values.any():
            df_by_lanes[df_name]["index"] = df_by_lanes[df_name]["index"].str.cat(df_by_lanes[df_name]["index2"], sep="-")
        df_by_lanes[df_name] = df_by_lanes[df_name].drop("index2", axis=1)
    # format two tables in the html with different column headers
    with open(write_to_file, 'w') as _file:
        _file.write("<h2>Lane Summary<h2>" + demux_stats_csv_convert.to_html(index = False, float_format =  '{:,.0f}'.format) + "\n<h2>Top Unknown Barcodes<h2>\n" + "<table>\n" )
        for value in df_by_lanes.values():
            _file.write("<td>" + value.to_html(index = False, float_format =  '{:,.0f}'.format) + "</td>")
        _file.write("\n</table>")

"""
demux_reports_dir - full path to the reports dir such as /igo/staging/FASTQ/DIANA_0442_BH3YJ7DSX3_V1/Reports/
run_folder - run folder only such as DIANA_0442_BH3YJ7DSX3
"""
def convert_dragen_reports(demux_reports_dir, run_folder):
    html_file = "/home/igo/html/{}_laneBarcode.html".format(run_folder)
    build_lane_summary_html(demux_reports_dir, html_file)
    #TODO set file creation time similar to:
    # touch $toNameLocal -r $dragen_replay #set correct timestamp on new html file from a DRAGEN demux file

if __name__ == '__main__':
    # Converting DRAGEN reports in folder /igo/staging/FASTQ/DIANA_0442_BH3YJ7DSX3/Reports/ to name /home/igo/html/DIANA_0442_BH3YJ7DSX3_laneBarcode.html
    #Usage: python dragen_csv_to_html.py [dragen_demux_dir] [output_file_name]
    build_lane_summary_html(sys.argv[1], sys.argv[2])