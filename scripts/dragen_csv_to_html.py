import pandas
import sys

if __name__ == "__main__":
    htmlfile = sys.argv[1]
    newhtmlstring = htmlfile.replace(".html", "_.html") # The name of the new HTML file to write to

    demux_stats_html = pandas.read_csv(htmlfile)

    demux_stats_html.to_html(newhtmlstring)