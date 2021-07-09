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
figureWidth=3

plt.figure(figsize=(figureWidth,figureHeight))

panel1Width=2
panel1Height=2
relativePanelWidth1=panel1Width/figureWidth
relativePanelHeight1=panel1Height/figureHeight

# main panel
panel1=plt.axes([0.1665,0.1665,relativePanelWidth1,relativePanelHeight1])
panel1.set_xlim(-12,12)
panel1.set_ylim(0,60)
plt.xlabel(r'$log_{2}$'"(fold change)")
plt.ylabel('-' + r'$log_{10}$' + '(p-value)')

# list for x, y values and name
all_fold=[]
all_pvals=[]
gene_names=[]

red_folds=[]
red_pvals=[]
label_folds=[]
label_pvals=[]
black_folds=[]
black_pvals=[]

# open specified input data file and parse information
infile=open(input_file,'r')

rand_list=[]
# scatter plot for main panel
for line in infile:
    try:
        splitLine=line.strip().split('\t')
        gene=splitLine[0]
        
        fold_val=(float)(splitLine[1])
        fold_abs=abs(fold_val)
        fold_change=(2**fold_abs)
        p_val=(float)(splitLine[2])
        y_val=((np.log10(p_val))*-1)
        # print("done try")
        if fold_change > 10 and y_val > 30 and fold_val < 0:
            label_folds.append(fold_val)
            label_pvals.append(y_val)
            gene_names.append(gene)
        elif fold_change > 10 and y_val > 8:
            red_folds.append(fold_val)
            red_pvals.append(y_val)
        else:
            black_folds.append(fold_val)
            black_pvals.append(y_val)
        
    except:
        continue

panel1.plot(red_folds,red_pvals, \
            marker='o',            
            markerfacecolor='red',          
            markeredgecolor='black',            
            markersize=1.45,            
            markeredgewidth=0,            
            linewidth=0)
      
panel1.plot(black_folds,black_pvals, \
            marker='o',            
            markerfacecolor='black',          
            markeredgecolor='black',            
            markersize=1.45,            
            markeredgewidth=0,            
            linewidth=0)

panel1.plot(label_folds,label_pvals, \
            marker='o',            
            markerfacecolor='red',          
            markeredgecolor='black',            
            markersize=1.45,            
            markeredgewidth=0,            
            linewidth=0)
for i, name in enumerate(gene_names):
    plt.annotate(name, (label_folds[i],label_pvals[i]), fontsize=6, xytext=(-1.42,-0.019), textcoords='offset points', ha='right', va='center')
    

plt.savefig(output_file,dpi=600)
