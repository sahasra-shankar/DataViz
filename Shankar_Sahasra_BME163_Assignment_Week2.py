import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument('-i', '--input_file')
parser.add_argument('-s', '--style_sheet', default='BME163')
parser.add_argument('-o', '--output_file')
args=parser.parse_args()
output_file=args.output_file
input_file=args.input_file

plt.style.use(args.style_sheet)

figureHeight=2
figureWidth=5

plt.figure(figsize=(figureWidth,figureHeight))

panel1Width=1
panel1Height=1
relativePanelWidth1=panel1Width/figureWidth
relativePanelHeight1=panel1Height/figureHeight

panel2Width=0.25
panel2Height=1
relativePanelWidth2=panel2Width/figureWidth
relativePanelHeight2=panel2Height/figureHeight

panel3Width=1
panel3Height=0.25
relativePanelWidth3=panel3Width/figureWidth
relativePanelHeight3=panel3Height/figureHeight

# main panel
panel1=plt.axes([0.17,0.15,relativePanelWidth1,relativePanelHeight1])
panel1.set_xlim(0,15)
panel1.set_ylim(0,15)
panel1.axes.yaxis.set_visible(False)

# side panel
panel2=plt.axes([0.1,0.15,relativePanelWidth2,relativePanelHeight2])
panel2.set_xlim(20,0)
panel2.set_ylim(0,15)

# top panel
panel3=plt.axes([0.17,0.69,relativePanelWidth3,relativePanelHeight3])
panel3.set_xlim(0,15)
panel3.axes.xaxis.set_visible(False)
panel3.set_ylim(0,20)

# list for x, y values and name
rep1_list=[]
rep2_list=[]
gene_names=[]

# open specified input data file and parse information
infile=open(input_file,'r')

for line in infile:
    splitLine=line.strip().split('\t')
    x_val=int(splitLine[1])
    y_val=int(splitLine[2])
    rep1_list.append(x_val)
    rep2_list.append(y_val)

log_rep1=np.log2(np.array(rep1_list)+1) #log array for rep1(x)
log_rep2=np.log2(np.array(rep2_list)+1) #log array for rep2(y)

# scatter plot for main panel
panel1.scatter(log_rep1,log_rep2, \
               facecolor=(0,0,0),\
               edgecolor='black',\
               s=2,\
               linewidth=0,\
               alpha=0.1)

bins= np.arange(0,15,0.5)
xHisto,bins2=np.histogram(log_rep1,bins)
yHisto,bins2=np.histogram(log_rep2,bins)

# histogram for top panel
for i in np.arange(0,len(xHisto),1):
    bottom=0
    left=bins[i]
    width=bins[i+1]-left
    height=np.log2(xHisto[i]+1)
    rectangle=mplpatches.Rectangle([left,bottom],width,height,                               
                                    facecolor=(0.5,0.5,0.5),                               
                                    edgecolor='black',                               
                                    linewidth=0.4)
    panel3.add_patch(rectangle)

#histogram for side panel
for i in np.arange(0,len(yHisto),1):
    bottom=0
    left=bins[i]
    width=bins[i+1]-left
    height=np.log2(yHisto[i]+1)
    rectangle=mplpatches.Rectangle([bottom,left],height,width,                               
                                    facecolor=(0.5,0.5,0.5),                               
                                    edgecolor='black',                               
                                    linewidth=0.4)
    panel2.add_patch(rectangle)

plt.savefig(output_file,dpi=600)


