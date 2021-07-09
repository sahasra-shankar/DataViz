import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib.image as img
import random
import sys
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument('-i', '--input_file')
parser.add_argument('-p', '--pngs')
parser.add_argument('-s', '--style_sheet', default='BME163')
parser.add_argument('-o', '--output_file')
args=parser.parse_args()
output_file=args.output_file
input_file=args.input_file
pngs_path=args.pngs

plt.style.use(args.style_sheet)

figureHeight=3
figureWidth=6

plt.figure(figsize=(figureWidth,figureHeight))

panel1Width=2.4
panel1Height=1
relativePanelWidth1=panel1Width/figureWidth
relativePanelHeight1=panel1Height/figureHeight

# left panel
panel1=plt.axes([0.10,0.3,relativePanelWidth1,relativePanelHeight1])
panel1.set_xlim(0,20)
panel1.set_ylim(0.0,2.0)
panel1.set_xticklabels(np.arange(-10, 11, 5))
panel1.set_xlabel('Distance to\nSplice Site')
panel1.set_ylabel('Bits')
panel1.set_title('5\'SS')

# right panel
panel2=plt.axes([0.57,0.3,relativePanelWidth1,relativePanelHeight1])
panel2.set_xlim(0,20)
panel2.set_ylim(0.0,2.0)
panel2.set_xticklabels(np.arange(-10, 11, 5))
panel2.axes.yaxis.set_visible(False)
panel2.set_xlabel('Distance to\nSplice Site')
panel2.set_ylabel('Bits')
panel2.set_title('3\'SS')

panel1.plot([10, 10], [0, 2], lw=1/2, color='black')
panel2.plot([10, 10], [0, 2], lw=1/2, color='black')

A = img.imread(pngs_path + "A_small.png")
T = img.imread(pngs_path + "T_small.png")
C = img.imread(pngs_path + "C_small.png")
G = img.imread(pngs_path + "G_small.png")

class FastAReader:
    def __init__ (self, fname=''):
        self.fname = fname
            
    def doOpen (self):
        if self.fname=='':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        header = ''
        sequence = ''
        with self.doOpen() as fileH:
            header = ''
            sequence = ''
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence
    
file_read=FastAReader(input_file)
file_read.doOpen()

headSeqDict={}
upList=[] #starts with 5
downList=[] #starts with 3

for head,seq in file_read.readFasta():
    headSeqDict[head]=seq
    if head.split("_")[0]=="3\'":
        downList.append(seq)
    else:
        upList.append(seq)

upSeqCount=len(upList)
downSeqCount=len(downList)

# Get base frequency based on base counts at each bp (0-20bp)
def base_freq(stream_list,seq_count):
    def base_counts(stream_list):
        BClist=[]
        for i in range(20):
            BClist.append({'A':0,'C':0,'G':0,'T':0})
        for seq in stream_list:
            for pos in range(20):
                if seq[pos]=='A':
                    BClist[pos]["A"] +=1
                elif seq[pos]=='C':
                    BClist[pos]["C"] +=1
                elif seq[pos]=='G':
                    BClist[pos]["G"] +=1
                elif seq[pos]=='T': 
                    BClist[pos]["T"] +=1
        return BClist

    BCdictList=base_counts(stream_list)
    BFlist=[]
    for i in range(20):
        BFlist.append({'A':0,'C':0,'G':0,'T':0})
    for countDictIndex in range(len(BCdictList)):
        for nuc,count in BCdictList[countDictIndex].items():
            BFlist[countDictIndex][nuc]=(count/seq_count)

    return BFlist

upBF=base_freq(upList,upSeqCount)
downBF=base_freq(downList,downSeqCount)

# Get base heights based on base frequency per bp
def base_heights(freq_dict_list,seq_count):
    all_heights=[]
    approxCorrection=(3/(np.log(2)*(2*seq_count)))
    for i in range(20):
        all_heights.append({'A':0,'C':0,'G':0,'T':0})
   
    for i in range(len(freq_dict_list)):
        total=0
        for base,freq in freq_dict_list[i].items():
            total+=((freq)*np.log2(freq))  
        uncertainty = (total*-1)
        for base,freq in freq_dict_list[i].items():
            infoContent=((np.log2(4))-(uncertainty+approxCorrection))
            baseHeight=freq*infoContent
            all_heights[i][base]=baseHeight

    return all_heights

upHeightDict=base_heights(upBF,upSeqCount)
downHeightDict=base_heights(downBF,downSeqCount)

# Sort base heights from smallest to largest and plot in respective panels
def sorted_heights(height_dict_list,panel):
    sHeightList=[]
    for hIndex in range(len(height_dict_list)):
        sorted_hdict = {}
        sorted_keys = sorted(height_dict_list[hIndex], key=height_dict_list[hIndex].get) 
        for val in sorted_keys:
            sorted_hdict[val] = height_dict_list[hIndex][val]
        sHeightList.append(sorted_hdict)
        
        start=0
        for key, height in sHeightList[hIndex].items():
            baseBot=start
            baseTop=baseBot+height
            start=baseTop

            if key == 'A':
                panel.imshow(A, aspect='auto', extent=[hIndex, hIndex+1,
                            baseBot, baseTop])
            elif key == 'G':
                panel.imshow(G, aspect='auto', extent=[hIndex, hIndex+1,
                            baseBot, baseTop])
            elif key == 'C':
                panel.imshow(C, aspect='auto', extent=[hIndex, hIndex+1,
                            baseBot, baseTop])
            elif key == 'T':
                panel.imshow(T, aspect='auto', extent=[hIndex, hIndex+1,
                            baseBot, baseTop])
    return 

sorted_heights(upHeightDict,panel1)
sorted_heights(downHeightDict,panel2)

plt.savefig(output_file,dpi=600)