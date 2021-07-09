import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import random
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
figureWidth=7

plt.figure(figsize=(figureWidth,figureHeight))

panel1Width=5
panel1Height=2
relativePanelWidth1=panel1Width/figureWidth
relativePanelHeight1=panel1Height/figureHeight

# main panel
panel1=plt.axes([0.1,0.2,relativePanelWidth1,relativePanelHeight1])
panel1.set_xlim(0.25,11.75)
panel1.set_ylim(75,100)
panel1.set_xlabel('Subread Coverage')
panel1.set_ylabel('Identity %')

locs = np.arange(1,12)
labels = ('1','2','3','4','5','6','7','8','9','10','>10')
plt.xticks(locs,labels)
panel1.plot([0,12], [95,95], linestyle=(0, (1, 1.9, 2, 1.9)), linewidth=0.5, color='black')
#plt.axhline(x=0,y=95,color='black',linestyle=(0, (3, 5, 1, 5)),linewidth=1.5)

#subcov_bins=[]
all_ids=[]

coverageDict={}
#initialize each index of dict to be list
for bins in range(1,12,1):
    coverageDict[bins]=[]

# open specified input data file and parse information
infile=open(input_file,'r')
for line in infile:
    splitLine=line.strip().split('\t')
    full_read_name=splitLine[0]
    cov_bin=(int)(full_read_name.split("_")[3])
    id=(float)(splitLine[1])
    if cov_bin < 11:
        coverageDict[cov_bin].append(id)
    else:
        coverageDict[11].append(id)
    
# Custom Swarm Function
def create_swarm(yList, xCenter, pointSize, xRange, yRange, panWidth, panHeight):
    pointsPlaced=[]
    minDist=pointSize/72

    switchCount=0 #also works
    for y_id in yList:
        subCenter=xCenter
        if(len(pointsPlaced)==0):
            pointsPlaced.append((subCenter,y_id))
        else:
            placed=False
            #switchCount=0 #kind of works
            for i in pointsPlaced:
                totalDist=np.sqrt(((((subCenter-i[0])/xRange)*panel1Width)**2)+((((y_id-i[1])/yRange)*panel1Height)**2))
                if (totalDist>minDist):
                    placed=True
                else:
                    placed=False
                if (placed==False):
                    #for i in np.arange(0,0.5,0.001):
                    if (switchCount%2==0):
                        subCenter+=(minDist*3) #xpos shift  
                    else:  
                        subCenter-=(minDist*3)   
                switchCount+=1 #works
            pointsPlaced.append((subCenter,y_id))

    return [val for val in pointsPlaced]         

for key,val in coverageDict.items():
    sample=random.sample(val,1000) #sample data
    finalPoints=create_swarm(sample,key,0.9,11.50,25,panel1Width,panel1Height)

    for point in finalPoints:
        panel1.plot(point[0], point[1],
                marker='o',
                markerfacecolor='black',
                markeredgecolor='black',
                markersize=0.85,
                markeredgewidth=0,
                linewidth=0)

# MEDIANS
all_medians = []
for scbin in range(1, 12, 1):
    all_medians.append(np.median(coverageDict[scbin]))

for i in range(0, len(all_medians), 1):
    mid= i+1
    width= 0.3
    bottom= all_medians[i]
    panel1.plot([mid-width, mid+width], [bottom, bottom], lw=1,
                color='red')


plt.savefig(output_file,dpi=600)
