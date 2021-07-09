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

figureHeight=3
figureWidth=5

plt.figure(figsize=(figureWidth,figureHeight))

panel1Width=0.75
panel1Height=2.5
relativePanelWidth1=panel1Width/figureWidth
relativePanelHeight1=panel1Height/figureHeight

# main panel
panel1=plt.axes([0.10,0.10,relativePanelWidth1,relativePanelHeight1])
panel1.set_xlim(-0.5,7.5)
panel1.set_xticks([0,1,2,3,4,5,6,7])
panel1.set_xticklabels(['0','','6','','12','','18',''])
panel1.set_ylim(0,1262)
panel1.set_xlabel('CT')
panel1.set_ylabel('Number of Genes')

geneList=[]
with open(input_file,'r') as infile:
    next(infile)
    for line in infile:
        split = line.strip().split('\t')
        geneList.append([int(split[4]),int(split[5]),int(split[6]),int(split[7]),int(split[8]),int(split[9]),
                          int(split[10]),
                          int(split[11]),
                          float(split[13])])

geneList=sorted(geneList,key=lambda glist:glist[-1], reverse=True)

# Colors
R=np.linspace(255/255, 56/255, 101)
G=np.linspace(225/255, 66/255, 101)
B=np.linspace(40/255, 157/255, 101)

data=[0,0,0,0,0,0,0,0]

Y=0
for item in geneList:
    X=-0.5
    for i in range(8):
        data[i]=item[i]
    data=np.array(data)
    modArray = ((data-min(data))/(max(data)-min(data)))*100
    for val in modArray:
        val=int(val)
        rectangle=mplpatches.Rectangle([X, Y],3,1,linewidth=0,facecolor=(R[val],G[val],B[val]))
        panel1.add_patch(rectangle)
        X+=1
    Y+=1

plt.savefig(output_file,dpi=600)