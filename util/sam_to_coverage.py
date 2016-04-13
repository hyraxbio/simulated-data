import pysam
import plotly.plotly as py
from plotly.graph_objs import *
import sys


def plot(coverage):
    py.sign_in("imogen", "mtf1tawct2")

    traces = []
    traces.append(
        Scatter(
            x=range(len(coverage)),
            y=coverage,
            name = 'Coverage',
            line = Line(
                color = 'blue',
                width = 2.0),    
        )
    )

    data = Data(traces)
    layout = Layout(
        title="Coverage over gene",
        titlefont = Font(
            size = 24),
        yaxis=YAxis(
            title='Coverage',
            showgrid=False,
            autorange = True,
            titlefont = Font(
                size = 20)
        ),
        xaxis = XAxis(
            title = "Nucleotide",
            titlefont = Font(
                size = 20)
        )
       
    )

    fig = Figure(data=data, layout=layout)
    py.plot(fig, width=1024, height=400)


def getcoverage(in_file, ref_length):

    sam_file = pysam.AlignmentFile(in_file, "r")
    

    coverage = [0] * ref_length

    for sam in sam_file.fetch():

        if sam.is_unmapped:
            continue

        rstart = sam.reference_start
        rend = sam.reference_end

        #print rstart, rend

        coverage[rstart:rend] = [c + 1 for c in coverage[rstart:rend]]

    plot(coverage)


if __name__ == "__main__":
  getcoverage(sys.argv[1], int(sys.argv[2]))





