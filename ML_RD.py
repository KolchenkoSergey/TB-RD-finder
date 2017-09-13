from optparse import OptionParser
import os
import sys
from distutils.spawn import find_executable
import subprocess
import pysam
import pandas as pd
import numpy as np
from itertools import islice
import math
import time
import datetime
from collections import Counter
import pickle
def predict_group(df):
 predicted_groups = []
 if df.loc[df['RD']=='RD105','pred'].item()==1 and \
 df.loc[df['RD']=='RD142','pred'].item()==0 and \
 df.loc[df['RD']=='RD150','pred'].item()==0 and \
 df.loc[df['RD']=='RD181','pred'].item()==0 and \
 df.loc[df['RD']=='RD207','pred'].item()==0:
  predicted_groups.append('2.1')
 elif df.loc[df['RD']=='RD105','pred'].item()==1 and \
 df.loc[df['RD']=='RD181','pred'].item()==0 and \
 df.loc[df['RD']=='RD207','pred'].item()==1 and \
 df.loc[df['RD']=='RD142','pred'].item()==0 and \
 df.loc[df['RD']=='RD150','pred'].item()==0:
  predicted_groups.append('2.2')
 elif df.loc[df['RD']=='RD105','pred'].item()==1 and \
 df.loc[df['RD']=='RD181','pred'].item()==1 and \
 df.loc[df['RD']=='RD207','pred'].item()==1 and \
 df.loc[df['RD']=='RD142','pred'].item()==0 and \
 df.loc[df['RD']=='RD150','pred'].item()==0:
  predicted_groups.append('2.2.1')
 elif df.loc[df['RD']=='RD105','pred'].item()==1 and \
 df.loc[df['RD']=='RD150','pred'].item()==1 and \
 df.loc[df['RD']=='RD181','pred'].item()==1 and \
 df.loc[df['RD']=='RD207','pred'].item()==1 and \
 df.loc[df['RD']=='RD142','pred'].item()==0:
  predicted_groups.append('2.2.1.1')
 elif df.loc[df['RD']=='RD105','pred'].item()==1 and \
 df.loc[df['RD']=='RD150','pred'].item()==0 and \
 df.loc[df['RD']=='RD181','pred'].item()==1 and \
 df.loc[df['RD']=='RD207','pred'].item()==1 and \
 df.loc[df['RD']=='RD142','pred'].item()==1:
  predicted_groups.append('2.2.1.2')
 elif df.loc[df['RD']=='RD105','pred'].item()==1 and \
 df.loc[df['RD']=='RD181','pred'].item()==0 and \
 df.loc[df['RD']=='RD207','pred'].item()==1 and \
 df.loc[df['RD']=='RD142','pred'].item()==0 and \
 df.loc[df['RD']=='RD150','pred'].item()==0: 
  predicted_groups.append('2.2.2')
 elif  df.loc[df['RD']=='RD239','pred'].item()==1:
  predicted_groups.append('1')
 elif df.loc[df['RD']=='RD750','pred'].item()==1:
  predicted_groups.append('3')
 elif df.loc[df['RD']=='RD183','pred'].item()==1:
  predicted_groups.append('4.1.1.1')
 elif df.loc[df['RD']=='RD193','pred'].item()==1:
  predicted_groups.append('4.1.1.3')
 elif df.loc[df['RD']=='RD182','pred'].item()==1:
  predicted_groups.append('4.1.2.1')
 elif df.loc[df['RD']=='RD761','pred'].item()==1:
  predicted_groups.append('4.3.2.1')
 elif df.loc[df['RD']=='RD115','pred'].item()==1:
  predicted_groups.append('4.3.3')
 elif df.loc[df['RD']=='RD174','pred'].item()==1:
  predicted_groups.append('4.3.4')
 elif df.loc[df['RD']=='RD122','pred'].item()==1:
  predicted_groups.append('4.5')
 elif df.loc[df['RD']=='RD724','pred'].item()==1:
  predicted_groups.append('4.6.1')
 elif df.loc[df['RD']=='RD726','pred'].item()==1:
  predicted_groups.append('4.6.2')
 elif df.loc[df['RD']=='RD219','pred'].item()==1:
  predicted_groups.append('4.8')
 elif df.loc[df['RD']=='RD711','pred'].item()==1:  
  predicted_groups.append('5')
 elif df.loc[df['RD']=='RD702','pred'].item()==1:
  predicted_groups.append('6')
 else:
        #print('bizare')
  predicted_groups.append('unclassified')
 return(predicted_groups)
def process_region(bamfile, region_start,region_end,average_coverage_whole):
    #start_time = time.time()
    bamFP = pysam.Samfile(bamfile, "rb");
    reads = bamFP.fetch(bamFP.references[0],(region_start-1),region_end)
    mapq = []
    flags1 = 0
    flags2 = 0
    isize = []
    soft_clipping_right=[]
    soft_clipping_left=[]
    edit_distance_total = []
    nreads = 0
    coverage = 0
    for i in range(4):
        coverage = coverage + sum(bamFP.count_coverage(bamFP.references[0],(region_start-1),region_end)[i])
    coverage = coverage / float(region_end-region_start+1)    
    for read in reads:
        nreads = nreads + 1
        mate_is_unmapped = 1 if read.mate_is_unmapped else 0
        improper_pair = 0 if read.is_proper_pair else 1
        flags1 = flags1  + mate_is_unmapped
        flags2 = flags2 + improper_pair
        isize.append(abs(read.template_length))
        if len(read.cigar)!=0:
            if read.cigar[0][0]==4:
                soft_clipping_left.append(read.cigar[0][1])
            if read.cigar[-1][0]==4:
                soft_clipping_right.append(read.cigar[-1][1])
        if read.has_tag('NM'):
            edit_distance_total.append(read.get_tag('NM'))
        mapq.append(read.mapping_quality)
    #print("--- %s seconds ---" % (time.time() - start_time))
    coverage = coverage / float(average_coverage_whole)
    if nreads==0:
        nreads = 1
    return([coverage,np.nanmean(mapq),np.nanmedian(isize),flags1/float(nreads),flags2/float(nreads),np.nanmean(edit_distance_total),np.nanmean(soft_clipping_left),np.nanmean(soft_clipping_right)])
def is_tool(name):
 return(find_executable(name) is not None)
parser = OptionParser()
parser.add_option("-f", "--file", action="store",type="string",
                  help="Specify filename", dest="filename")
parser.add_option("-i","--input",action="store",type="string",help="Specify input folder",
                dest="input_folder",default=str(os.getcwd()))
parser.add_option("-m","--mode",dest="mode"
                  ,help="interaction mode: if bam, process starts directly from analyzing bam file, bam file should be sorted and indexed; if fastq, files will be processed with bowtie2 first, will look for files _1.fastq.gz and _2.fastq.fz; for fastq --refi and --reff should be specified")
parser.add_option("-o","--output",dest="output",type="string",action="store",help="outuput folder",default=str(os.getcwd()))
parser.add_option("--refi",dest="refi",help="Path to folder and basename",action="store",default="ref/NC_000962.3")
parser.add_option("--reff",dest="reff",help="Path to reference fasta",action="store",default="ref/NC_000962.3.fasta")
# TO DO: add custom RD list
(option,args) = parser.parse_args()
print("Speicifed filename: "+option.filename)
print("Specified folder: "+option.input_folder)
print("Specified mode: "+option.mode)
if option.mode=='fastq':
 print("Will look for files "+option.input_folder+"/"+option.filename+"_1.fastq.gz and "
         +option.input_folder+"/"+option.filename+"_2.fastq.gz")
 fastq1_file = option.input_folder+"/"+option.filename+"_1.fastq.gz"
 fastq2_file = option.input_folder+"/"+option.filename+"_2.fastq.gz"
 if os.path.exists(fastq1_file):
  print("Fastq 1 exists")
  fastq1_file_exist = True
 else:
  print("Fastq1 does not exist")
  fastq1_file_exist = False
 if os.path.exists(fastq1_file):
  print("Fastq 2 exists")
  fastq2_file_exist = True
 else:
  print("Fastq 2 does not exist")
  fastq2_file_exist = False
 if is_tool('bowtie2'):
  print("Bowtie2 found")
  bowtie2_found = True
 else:
  print("Bowtie2 not found")
  bowtie_found= False
 if is_tool('samtools'):
  print('Samtools found')
  samtools_found = True
 else:
  print("samtools not found")
  samtools_found = False
 if option.refi == None:
  print("refi not specified")
  refi_specified = False
 else: 
  refi_specified = True
 if option.reff == None:
  print('reff not specified')
  reff_specified = False
 else:
  reff_specified = True
 if not reff_specified or not refi_specified or not bowtie2_found or not samtools_found:
  sys.exit(1)
 subprocess.call(['bowtie2', '-q','-t','--local','-x',option.refi,'-1',fastq1_file,'-2',fastq2_file,'-S',option.output+'/'+option.filename+'.sam','--threads','6'])
 subprocess.call(['samtools','view','-bS' '-o',option.output+'/'+option.filename+".bam",option.output+'/'+option.filename+'.sam'])
 subprocess.call(['samtools', 'sort',option.output+'/'+option.filename+'.bam',option.output+'/'+option.filename])
 subprocess.call(['samtools', 'index',option.output+'/'+option.filename+'.bam'])
 os.remove(option.output+'/'+option.filename+'.sam')
 bam_file = option.output+'/'+option.filename+'.bam'
 bai_file = option.output+'/'+option.filename+'.bam.bai'
 files_OK = True

if option.mode=='bam':
 print('Will look for files '+option.input_folder+"/"+option.filename+".bam and "
         +option.input_folder+"/"+option.filename+".bam.bai")
 bam_file = option.input_folder+"/"+option.filename+".bam"
 bai_file = option.input_folder+"/"+option.filename+".bam.bai"
 if os.path.exists(bam_file):
  print("Bam file exists")
  bam_file_exist = True
 else:
  print("Bam file does not exist")
  bam_file_exist = False
 if os.path.exists(bai_file):
  print("Bai file exists")
  bai_file_exist = True
 else:
  print("Bai file does not exist")
  bai_file_exist = False	
 if not bai_file_exist or not bam_file_exist:
  print('File not found!')
  sys.exit(1)
 print('Files OK')
 files_OK = True
# TO DO: check if sorted
if files_OK:
#read RDs list
 RD_data = pd.read_csv('RD_data.csv',sep='\t')
 left_RD = RD_data['Start']
 right_RD = RD_data['Stop']
 RD_list = RD_data['RD']
 try:
  file = pysam.AlignmentFile(bam_file,'rb')
 except ValueError:
  print('Bad file')
  file_correct = False
  sys.exit(1)
 else:
  output = subprocess.check_output("samtools depth "+str(bam_file)+" |  awk '{sum+=$3} END { print sum/NR}'", shell=True)
  coverage = float(output.strip('\n'))
  print(coverage)
  file_correct = True
 if file_correct:
  df_total = pd.DataFrame()
  for j in range(len(RD_list)):
   df_RD = pd.DataFrame(data = [process_region(bam_file,int(left_RD[j]),int(right_RD[j]),coverage)])
   df_RD.columns = ['cov','mapq','isize','flags1','flags2','edistance','softclipleft','softclipright']
   df_RD['Sample'] = option.filename
   df_RD['RD'] = RD_list[j]
   df_total = pd.concat([df_total,df_RD])
  with open('model_cat.pickle') as f:
   model_cat = pickle.load(f)
  f.close()
  X = df_total.drop(['isize','Sample','RD','flags1','flags2'],axis=1)
  X.fillna(0,inplace=True)
  y = model_cat.predict(X)
  df_total['pred'] = y
  df_total = df_total.drop_duplicates(subset=['RD','Sample'])
  pred_groups = predict_group(df_total)
  print(pred_groups)
  df_total.to_csv(option.filename+'.csv',sep='\t')
  f = open(option.filename+'.group','w')
  f.write(option.filename+'\t'+' '.join(pred_groups)+'\n')
  f.close() 
