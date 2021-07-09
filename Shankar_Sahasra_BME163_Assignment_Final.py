import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('-i1', '--input_file1')
parser.add_argument('-i2', '--input_file2')
parser.add_argument('-g', '--gtf_file')
parser.add_argument('-o', '--output_file',default='Shankar_Sahasra_BME162_Assignment_Final.png')
args=parser.parse_args()
output_file=args.output_file
input_file1=args.input_file1
input_file2=args.input_file2
gtf_file=args.gtf_file

plt.style.use('BME163')

figureHeight=5
figureWidth=10

plt.figure(figsize=(figureWidth,figureHeight))

panelWidth=10
panelHeight=1.25

relativePanelWidth=panelWidth/figureWidth
relativePanelHeight=panelHeight/figureHeight

# B->Input_5.psl; M->Input_6.psl; T->gencode.vM12.annotation.gtf
panelB=plt.axes([0,0.05,relativePanelWidth,relativePanelHeight])
panelM=plt.axes([0,0.35,relativePanelWidth,relativePanelHeight])
panelT=plt.axes([0,0.65,relativePanelWidth,relativePanelHeight])

# Coordinates for gene locus: chromosome(chr7), left = 45232945, right = 45240000
target=['chr7',45232945,45240000]

def readPSL(pslfile):
    readList=[]
    oFile=open(pslfile, 'r')
    for line in oFile:
        a=line.strip().split('\t')
        chromosome=a[13]
        start=int(a[15])
        end=int(a[16])
        blockstarts=np.array(a[20].split(',')[:-1],dtype=int)
        blockwidths=np.array(a[18].split(',')[:-1],dtype=int)
        read=[chromosome,start,end,blockstarts,blockwidths,False]
        readList.append(read)
    return readList


def plotReads(panel, readList, target):
    genome_chrom,genome_start,genome_end=target[0],target[1],target[2]
    readsToPlotList=[]
    for read in readList:
        chromosome,start,end,blockstarts,blockwidths,plotted=read[0],read[1],read[2],read[3],read[4],read[5]

        if chromosome==genome_chrom:
            if genome_start<start<genome_end or genome_start<end<genome_end:
                readsToPlotList.append(read)

    for yPos in range(1,len(readsToPlotList),1):
        end_of_read=0
        for read in readsToPlotList:
            if read[5] is False:
                bottom=yPos
                chromosome,start,end,blockstarts,blockwidths,plotted=read[0],read[1],read[2],read[3],read[4],read[5]
                if start>end_of_read:
                    rectangle=mplpatches.Rectangle([start,yPos+0.18],end-start,0.1,
                                                     facecolor='black',
                                                     edgecolor='black',
                                                     linewidth=0)
                    panel.add_patch(rectangle)

                    for index in np.arange(0,len(blockstarts),1):
                        blockstart=blockstarts[index]
                        blockwidth=blockwidths[index]
                        if len(read)==7:
                            if read[6][index]=='exon':
                                rectangle=mplpatches.Rectangle([blockstart,yPos+0.1],blockwidth,0.25,
                                                                 facecolor='black',
                                                                 edgecolor='black',
                                                                 linewidth=0)
                            else:
                                rectangle = mplpatches.Rectangle([blockstart,yPos],blockwidth,0.5,
                                                                 facecolor='black',
                                                                 edgecolor='black',
                                                                 linewidth=0)
                        else:
                            rectangle=mplpatches.Rectangle([blockstart,yPos],blockwidth,0.5,
                                                                facecolor='black',
                                                                edgecolor='black',
                                                                linewidth=0)
                        panel.add_patch(rectangle)
                    end_of_read=end
                    read[5]=True

    return bottom

transcriptList=[]
gtfDict={}
file=open(gtf_file,'r')
lines=file.readlines()[5:]
for line in lines:
    a=line.strip().split('\t')
    chromosome=a[0]
    read_type=a[2]

    if read_type=='exon' or read_type=='CDS':
        start=int(a[3])
        end=int(a[4])
        transcript=a[8].split(' transcript_id "')[1].split('"')[0]
        
        if transcript not in gtfDict.keys():
            gtfDict[transcript]=[]
        gtfDict[transcript].append([chromosome,start,end,read_type])
    #print(gtfdict)

for transcript,info in gtfDict.items():
    starts=[]
    ends=[]
    blockstarts=[]
    blockwidths=[]
    types=[]
    for detail in info:
        starts.append(detail[1])
        ends.append(detail[2])
        blockstarts.append(detail[1])
        blockwidths.append(detail[2]-detail[1])
        chromosome=detail[0]
        types.append(detail[3])
    transcriptList.append([chromosome,min(starts),max(ends),blockstarts,blockwidths,False,types])


sorted_reads1=sorted(readPSL(input_file1),key=lambda a:a[2])
sorted_reads2=sorted(readPSL(input_file2),key=lambda a:a[2])
sorted_readList3=sorted(transcriptList,key=lambda a:a[2])


bottom1=plotReads(panelB,sorted_reads1,target)
bottom1=bottom1+(bottom1*0.1)
panelB.set_xlim(target[1],target[2])
panelB.set_ylim(0,bottom1)
panelB.axes.xaxis.set_visible(False)
panelB.axes.yaxis.set_visible(False)

bottom2=plotReads(panelM,sorted_reads2,target)
bottom2=bottom2+(bottom2*0.1)
panelM.set_xlim(target[1],target[2])
panelM.set_ylim(0,bottom2)
panelM.axes.xaxis.set_visible(False)
panelM.axes.yaxis.set_visible(False)

bottom3=plotReads(panelT, sorted_readList3, target)
panelT.set_ylim(0,2+bottom3)
panelT.set_xlim(target[1],target[2])
panelT.axes.xaxis.set_visible(False)
panelT.axes.yaxis.set_visible(False)

plt.savefig(output_file, dpi=1200)
