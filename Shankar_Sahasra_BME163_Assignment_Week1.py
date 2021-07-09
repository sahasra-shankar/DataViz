import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument('-s', '--style_sheet', default='BME163')
parser.add_argument('-o', '--output_file')
args=parser.parse_args()
output_file=args.output_file

plt.style.use(args.style_sheet)

figureHeight=2
figureWidth=3.42

plt.figure(figsize=(figureWidth,figureHeight))

panelWidth=1
panelHeight=1
relativePanelWidth=panelWidth/figureWidth
relativePanelHeight=panelHeight/figureHeight

def apply_func(valArr):
     xList = []
     yList = []
     for val in valArr:
         if np.cos(val) > -0.05 and np.cos(val) < 1.05 and np.sin(val) > -0.05 and np.sin(val) < 1.06:
            xList.append(np.cos(val))
            yList.append(np.sin(val))
     return xList, yList

### PANEL 1
# left,bottom, width,height
panel1=plt.axes([0.13,0.22,relativePanelWidth,relativePanelHeight])
panel1.axes.xaxis.set_visible(False)
panel1.axes.yaxis.set_visible(False)

black=(0,0,0)
white=(1,1,1)

R=np.linspace(black[0],white[0],23)
G=np.linspace(black[1],white[1],23)
B=np.linspace(black[2],white[2],23)

xList=np.arange(0,(np.pi)/2,0.065)
xValList, yValList = apply_func(xList)

for c in range(0,23,1):
    panel1.plot(yValList[c],xValList[c],            
        marker='o',            
        markerfacecolor=(R[c],G[c],B[c]),          
        markeredgecolor='black',            
        markersize=2,            
        markeredgewidth=0,            
        linewidth=0)
        #print(R[c]) 
panel1.set_xlim(-0.03,1.05)
panel1.set_ylim(0,1.05)

### PANEL 2
panel2=plt.axes([0.56,0.22,relativePanelWidth,relativePanelHeight])
panel2.axes.xaxis.set_visible(False)
panel2.axes.yaxis.set_visible(False)
for i in np.arange(0,1,0.1):
    for j in np.arange(0,1,0.1):
        rectangle=mplpatches.Rectangle([i,j],0.1,0.1,                               
                                    facecolor=(i,j,1),                               
                                    edgecolor='black',                               
                                    linewidth=1)
        panel2.add_patch(rectangle)

plt.savefig(args.output_file,dpi=600)
