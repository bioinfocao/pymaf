"""
util module for extraction of Fasta sequence from whole genome multiple alignments MAF file
"""
# Authors: Huojun Cao <bioinfocao at gmail.com>
# License: BSD 3 clause
import glob
import os.path
from Bio import SeqIO
from datetime import datetime
import subprocess
import re
ensm=re.compile(r"EN[a-zA-Z0-9_]+\.[0-9]{1,3}")
ncbi=re.compile(r"[NX][MPR]_[0-9]{3,20}\.[0-9]{1,3}")
oneTwoDots=re.compile(r"[a-zA-Z]+[a-zA-Z0-9_]+\.[0-9]{1,3}[\.0-9]{0,3}")
transcriptID_regObj_list=[ensm,ncbi,oneTwoDots]

mm10_multiz60way_29mammal_species=["mm10", "rn5", "speTri2", "hg19", "dipOrd1", "panTro4", \
                                   "rheMac3", "canFam3", "myoLuc2", "vicPac1", "equCab2", "pteVam1", \
                                   "tupBel1", "dasNov3", "micMur1", "otoGar3", "loxAfr3", "felCat5", \
                                   "tarSyr1", "oryCun2", "turTru2", "sorAra1", "bosTau7", "cavPor3", \
                                   "ochPri2", "choHof1", "echTel1", "proCap1", "eriEur1"]

hg38_multiz100way_58mammal_species=["hg38","panTro4","gorGor3","ponAbe2","nomLeu3","rheMac3","macFas5",\
                                    "papAnu2","chlSab2","calJac3","saiBol1","otoGar3","tupChi1","speTri2",\
                                    "jacJac1","micOch1","criGri1","mesAur1","mm10","rn6","hetGla2","cavPor3",\
                                    "chiLan1","octDeg1","oryCun2","ochPri3","susScr3","vicPac2","camFer1",\
                                    "turTru2","orcOrc1","panHod1","bosTau8","oviAri3","capHir1","equCab2",\
                                    "cerSim1","felCat8","canFam3","musFur1","ailMel1","odoRosDiv1","lepWed1",\
                                    "pteAle1","pteVam1","eptFus1","myoDav1","myoLuc2","eriEur2","sorAra2",\
                                    "conCri1","loxAfr3","eleEdw1","triMan1","chrAsi1","echTel2","oryAfe1","dasNov3"]

def getCmdArgs(argv):
    """
    Parse command line argvs
    Example: python test.py test.txt --useZip -b bfile
    pos_arg_list
    ['test.py', 'test.txt']
    named_arg_dict
    {'--useZip': None, '-b': 'bfile'}
    """
    pos_arg_list=[] # position arguments
    named_arg_dict = {}  # Empty dictionary to store key-value pairs.
    argv=[el.strip(" ") for el in argv]
    while argv:  # While there are arguments left to parse...
        if argv[0][0] == '-':  
            if argv[0][0:2] == '--': # Found a "--opt [value]"
                if len(argv)>1 and (argv[1][0] != "-"):
                    named_arg_dict[argv[0]] = argv[1]  # Add key and value to the dictionary.
                    argv = argv[2:] # Reduce the argument list by copying it starting from index 2.
                else:
                    named_arg_dict[argv[0]] = ""  # Add key and value to the dictionary.
                    argv = argv[1:]
            else: # Found a "-name value" pair.
                if len(argv)>1 and argv[1][0] != "-":
                    named_arg_dict[argv[0]] = argv[1]  # Add key and value to the dictionary.
                    argv = argv[2:] # Reduce the argument list by copying it starting from index 2.
                else:
                    named_arg_dict[argv[0]] = None
                    argv = argv[1:]
        else:
            pos_arg_list.append(argv[0])
            argv = argv[1:]  
    return pos_arg_list,named_arg_dict

def getArgsVal(named_arg_dict, name_list):
    val_list=[named_arg_dict.get(name,None) for name in name_list]
    return val_list

################################################################

def excu(cmdText, retErr=False, verbose=True):
    """ Call an external cmd and return result """
    if verbose:
        print("{0} Start:\n {1}".format(datetime.now(), cmdText))
    p = subprocess.Popen(cmdText,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out, err = p.communicate()
    if verbose:
        print("{0} Finished:\n {1}".format(datetime.now(), cmdText))
    if retErr:
        return out,err
    else:
        return out

################################################################

def getFiles(input_dir,fileNameEnd=".fasta"):
    """
    Get list of files
    """
    if len(fileNameEnd)>0 and fileNameEnd[0]!=".":fileNameEnd="."+fileNameEnd
    files=glob.glob(os.path.join(input_dir,"*"+fileNameEnd))
    return files

def getGeneName(infile,fileNameEnd):
    geneName=os.path.split(infile)[1].rstrip(fileNameEnd)
    return geneName

def colLineSPDict(infile,colIndex):
    """ for a given input file, return a dict, key is the selected col, value is linesp """
    col_lineSP_dict={}
    with open(infile,'r') as f:
        for line in f:
            if line.startswith("#"): continue
            linesp=line.rstrip().split("\t")
            if len(linesp[0])>5: continue # only keep standard chr
            col_lineSP_dict[linesp[colIndex]]=linesp
    return col_lineSP_dict

def getLinespList(infile,sep="\t"):
    linesp_list=[]
    with open(infile,'r') as f:
        for line in f:
            linesp=line.split(sep)
            linesp[-1]=linesp[-1].rstrip()
            linesp_list.append(linesp)
    return linesp_list

def writeLinespList(outfile,linespList,sep="\t"):
    txt="\n".join([sep.join([str(el) for el in linesp]) for linesp in linespList])
    with open(outfile,'w') as f:
        f.write(txt)

################################################################

def listToCount(aList,dsc=True,minCount=1):
    """ given a list return a list of sorted tuple (item,count) """
    countDict={}
    for item in aList:
        countDict[item]=countDict.get(item,0)+1
    if minCount>1:
        countList=sorted([(key,val) for key,val in countDict.items() if val>=minCount], key=lambda x:x[1], reverse=dsc)
    else:
        countList=sorted([(key,val) for key,val in countDict.items()], key=lambda x:x[1], reverse=dsc)
    return countList
        
################################################################################

def RNAToDNAPos(rnaPos, bed12Linesp):
    """ Given bed12 blockSizes, blockStarts, convert RNA pos to DNA pos """
    blockSizes=[int(el) for el in bed12Linesp[10].rstrip(",").split(",")]
    blockStarts=[int(el) for el in bed12Linesp[11].rstrip(",").split(",")]
    rnaLen=0
    dnaPos,relativeDNAPos=None,None
    for bStart,bSize in zip(blockStarts,blockSizes):
        if rnaPos <= rnaLen + bSize:
            relativeDNAPos=bStart + (rnaPos-rnaLen)
            break
        rnaLen+=bSize
    if relativeDNAPos is None:
        print("bed12Linesp:{0}\trnaPos:{1}".format(bed12Linesp,rnaPos))
    else:
        dnaPos=int(bed12Linesp[1]) + relativeDNAPos
    return dnaPos

def DNAtoRNAPos(dnaPos, bed12Linesp):
    """ Given bed12 blockSizes, blockStarts, convert DNA pos to RNA pos """
    blockSizes=[int(el) for el in bed12Linesp[10].rstrip(",").split(",")]
    blockStarts=[int(el) for el in bed12Linesp[11].rstrip(",").split(",")]
    rnaLen=0
    rnaPos=None
    relativeDNAPos=dnaPos-int(bed12Linesp[1])
    for bStart,bSize in zip(blockStarts,blockSizes):
        if relativeDNAPos <= bStart + bSize:
            rnaPos=rnaLen+(relativeDNAPos-bStart)
            break
        rnaLen+=bSize
    return rnaPos

def negativeToPositive(pos,totalLen):
    """
    Convert a position of negative strand to positive strand
    Use BED zero based format, end exclusive
    """
    positive_pos=totalLen-pos
    return positive_pos

def trimCodonSeq(fastaSeq,frame=0):
    """ Trim seq to make its length is 3X """
    outframe=len(fastaSeq[frame:]) % 3
    if outframe != 0:
        trimmedSeq=fastaSeq[frame:-outframe]
    else:
        trimmedSeq=fastaSeq[frame:]
    return trimmedSeq

def getFirstFastaSeq(fasta_file, ungap=True):
    with open(fasta_file,'r') as in_f:
        fasta_sequences = SeqIO.parse(in_f,'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, fasta.seq
            if ungap: sequence=sequence.ungap("-") # remove -
            return sequence
        
def getFirstLine(infile):
    with open(infile,'r') as f:
        line=f.readline()
    return line

################################################################################

def regTranscriptID(text):
    """
    find transcript id in the text by regex expression
    Some example:
    ENSMUST00000187631.1
    NM_025300.4
    MSTRG.6.5
    """
    for regObj in transcriptID_regObj_list:
        match_list=regObj.findall(text)
        if match_list:
            return match_list[0]

def getTranscriptID(text):
    """
    Find first transcript id in the text. Some example:
    transcript_id "NM_025300.4"
    gene_id "MSTRG.6.5"
    Name=CG11023
    """
    if text.find("transcript_id")>=0:
        transcript_id=text.split("transcript_id")[1].split(";")[0].split("-")[0].replace('"','').replace("'","").strip("=").strip()
    elif text.find("gene_id")>=0:
        transcript_id=text.split("gene_id")[1].split(";")[0].replace('"','').replace("'","").strip()
    elif text.find("Name=")>=0:
        transcript_id=text.split("Name=")[1].split(";")[0].strip()
    elif text.find("transcript:")>=0:
        transcript_id=text.split("transcript:")[1].split("-")[0].strip()
    else:
        transcript_id=regTranscriptID(text)
    return transcript_id

def getTranscriptIDSet(gtfFile):
    """
    extract all TranscriptIDs from gtfFile
    """
    linesp_list=getLinespList(gtfFile,sep="\t")
    transcriptID_list=[]
    for linesp in linesp_list:
        if len(linesp)==9:
            transcriptID=getTranscriptID(linesp[8])
            if transcriptID: transcriptID_list.append(transcriptID)
    transcriptID_set=set(transcriptID_list)
    return transcriptID_set

################################################################################

def gtfToBed6(infile,outfile,verbLevel=1):
    """
    Convert transcripts' GTF format to BED (6 columns) format
    Only keep exons. Seperate chr, and within each chr sorted based on start
    Basically take these columns from GTF file
    $0,$-1,$4,transcript_id,$5,$6
    For BED, It's the same regardless of the strand. The first base in the interval is at the start value and the last is end-1. The first base in the genome is 0.
    For GTF, it's the same regardless of strand, the first base in the interval is at the start value and the last is end. The first base in the genome is 1.
    GTF:start,end = BED:start-1,end
    """
    bed_line_list_dict={}
    if verbLevel>0:
        print("Convert GTF file:\n{0}\n to sorted BED6 file:\n{1}".format(infile,outfile))
    with open(infile,'r') as in_f:
        for line in in_f:
            linesp=line.rstrip().split("\t")
            # only keep lines that have 9 columns and exon in 3rd column
            if len(linesp)==9 and linesp[2].lower().find("exon")>=0:
                transcript_id=getTranscriptID(linesp[8])
                if not transcript_id: transcript_id=linesp[8]
                chrN=linesp[0]
                if chrN not in bed_line_list_dict: bed_line_list_dict[chrN]=[]
                bed_line=(chrN,int(linesp[3])-1,linesp[4],transcript_id+"_exon_",linesp[5],linesp[6])
                bed_line_list_dict[chrN].append(bed_line)
    # seperate chr, and within each chr sorted based on start
    all_sorted_bed_line_list=[]
    for chrN in bed_line_list_dict:
        bed_line_list=bed_line_list_dict[chrN]
        sorted_bed_line_list=sorted(bed_line_list,key=lambda x: x[1])
        all_sorted_bed_line_list+=sorted_bed_line_list
    # write to BED6
    all_sorted_bed_lineTexts=["\t".join([str(el) for el in line]) for line in all_sorted_bed_line_list]
    with open(outfile,'w') as out_f:
        out_f.write("\n".join(all_sorted_bed_lineTexts))
        
        
##############################################################################
## concatenate MAF of exons for each transcript
def getTranscriptExonsMAF(exon_bed_file,chr_start_end_maf_dir,filterChrNNN=True):
    """
    limit to standard chromosomes 1-22 + gender
    skip chrUn, random, alts
    """
    transcript_chrStartEndList_dict={}
    with open(exon_bed_file,'r') as f:
        for line in f:
            linesp=line.rstrip().split("\t")
            # limit to standard chromosomes + gender
            # chrN,chrNN,chrNNN are OK
            # chrNNN_xxx is removed  
            if filterChrNNN and len(linesp[0])>6: continue 
            if linesp[3].find("_exon_")>=0:
                transcriptName=linesp[3].split("_exon_")[0]
            else:
                transcriptName=getTranscriptID(linesp[3])
                if not transcriptName: transcriptName=linesp[3]
            chrN=linesp[0]
            if transcriptName not in transcript_chrStartEndList_dict: transcript_chrStartEndList_dict[transcriptName]=[]
            transcript_chrStartEndList_dict[transcriptName].append("_".join(linesp[:3]))
    transcript_mafFiles_list=[]
    for transcript,chrStartEndList in transcript_chrStartEndList_dict.items():
        sortedChrStartEndList=sorted(chrStartEndList, key=lambda chrStartEnd: chrStartEnd.split("_")[1])
        maf_files=[os.path.join(chr_start_end_maf_dir,chrStartEnd+".maf") for chrStartEnd in sortedChrStartEndList]
        transcript_mafFiles_list.append((transcript,maf_files))
        
    return transcript_mafFiles_list

def catExonsMAFCmds(exon_bed_file,chr_start_end_maf_dir,transcript_catMAF_dir,freshTransCatMAF,filterChrNNN):
    """
    cat 1.maf 2.maf > transcriptID.maf
    """
    transcript_mafFiles_list=getTranscriptExonsMAF(exon_bed_file,chr_start_end_maf_dir,filterChrNNN)
    cmd_list=[]
    for transcript, mafFiles in transcript_mafFiles_list:
        if not os.path.exists(transcript_catMAF_dir): os.makedirs(transcript_catMAF_dir) 
        transcriptMaf_file=os.path.join(transcript_catMAF_dir,transcript+".maf")
        if (not freshTransCatMAF) and os.path.isfile(transcriptMaf_file) and os.path.getsize(transcriptMaf_file) > 1: continue
        cmdText="cat {0} > {1}".format(" ".join(mafFiles),transcriptMaf_file)
        cmd_list.append(cmdText)
    return cmd_list

##############################################################################
## Modify and split exon bed6 file

def filterExistedMafRegions(linesp_list, maf_dir, smallestMafSize=0):
    filtered_linesp_list=[]
    for linesp in linesp_list:
        maf_file=os.path.join(maf_dir,"{0}_{1}_{2}.maf".format(linesp[0],linesp[1],linesp[2]))
        if os.path.exists(maf_file) and os.path.getsize(maf_file)>smallestMafSize:
            continue
        else:
            filtered_linesp_list.append(linesp)
    return filtered_linesp_list

def splitIntoNonoverlappingList(linesp_list):
    """
    Given a list of lines, split them into few list that in each list there is no overlapping regions.
    That is sort by start and the start is always bigger than last line's end
    """
    if len(linesp_list)==0: return []
    nonover_linespList_list=[]
    firstLine=linesp_list[0]
    # add first list to first linespList
    nonover_linespList_list.append([firstLine,])
    if len(linesp_list)<2: return nonover_linespList_list
    for linesp in linesp_list[1:]:
        inserted=False
        for linespList in nonover_linespList_list:
            last_end=linespList[-1][2]
            if linesp[1]>last_end:
                linespList.append(linesp)
                inserted=True
                break
        if not inserted:
            nonover_linespList_list.append([linesp,])
    return nonover_linespList_list
        

def exonBedSplitAndChangeName(infile,maf_dir,outfile_dir=""):
    """
    Modify and split exon bed6 file
    change name to chr_start_end
    Split exon bed based on chr and strand. Also split exon bed into
    non overlapping files
    So that mafsInRegion can works properly
    chr3	129213931	129214781	NM_001042502.2_exon_0_0_chr3_129213932_f	0	+
    To
    chr3	129213931	129214781	chr3_129213931_129214781	0	+
    Also check if the maf already existed chr_start_end.maf (chrY_90667524_90667625.maf)
    """
    chrStrand_linespList_dict={}
    with open(infile,'r') as inf:
        for line in inf:
            linesp=line.rstrip().split("\t")
            if len(linesp)==6:
                chrStrand="{0}{1}".format(linesp[0],linesp[5])
                if chrStrand not in chrStrand_linespList_dict: chrStrand_linespList_dict[chrStrand]=set()
                # change name to chr_start_end
                linesp[3]="{0}_{1}_{2}".format(linesp[0],linesp[1],linesp[2])
                # convert start,end to int
                linesp[1]=int(linesp[1])
                linesp[2]=int(linesp[2])                
                chrStrand_linespList_dict[chrStrand].add(tuple(linesp))
    #
    for chrStrand in chrStrand_linespList_dict:
        linespList=list(chrStrand_linespList_dict[chrStrand])
        sortedLinespList=sorted(linespList, key=lambda linesp: linesp[1])
        # remove regions that maf already extracted (maf file chr_start_end.maf already existed)
        # if maf_dir not exist or not defined, skip this step
        if os.path.isdir(maf_dir):
            sortedLinespList=filterExistedMafRegions(sortedLinespList, maf_dir)        
        # split lines into non overlapping list
        nonover_linespList_list=splitIntoNonoverlappingList(sortedLinespList)
        if not outfile_dir:
            outfile_dir=os.path.dirname(infile)        
        for index,linespList in enumerate(nonover_linespList_list):
            outfile=os.path.join(outfile_dir,chrStrand.replace("+","plus").replace("-","minus")+"_"+str(index)+".bed")
            linespList_text=["\t".join([str(el) for el in linesp]) for linesp in linespList ]
            with open(outfile,'w') as outf:
                outf.write("\n".join(linespList_text))

##############################################################################
## generate cmds for extraction of exons' MAF by mafsInRegion
                
def mafExtractCmds(mafsInRegion_bin,bed_files_dir,maf_dir,extracted_chr_start_end_maf_dir, gzInput=True ,verbLevel=1):
    """
    [gunzip -c chr3.maf.gz > chr3.maf]
    mafsInRegion noncoding_exons_mm10_VM16_chr3minus_2.bed -outDir mm10_multiz60way_chr_start_end_maf chr3.maf
    [rm chr3.maf]
    """
    splitted_bed_files=glob.glob(os.path.join(bed_files_dir,"*.bed"))
    chrN_cmd_list={}
    for bed_file in splitted_bed_files:
        headPath, tailFN = os.path.split(bed_file)
        if tailFN.find("plus")>=0:
            chrN=tailFN.split("plus")[0]
        else:
            chrN=tailFN.split("minus")[0]
        cmd="{0} {1} -outDir {2} {3}".format(mafsInRegion_bin, bed_file, extracted_chr_start_end_maf_dir, os.path.join(maf_dir,chrN+".maf"))
        if verbLevel>0:
            cmd="echo 'Working on: {0}'\n".format(bed_file) +cmd
        if chrN not in chrN_cmd_list: chrN_cmd_list[chrN]=[]
        chrN_cmd_list[chrN].append(cmd)
    # combine and add gunzip support
    combined_cmd_list=[]
    for chrN in chrN_cmd_list:
        cmd_list=chrN_cmd_list[chrN]
        if gzInput:
            gunzip_cmd="gunzip -c {0} > {1}".format(os.path.join(maf_dir,chrN+".maf.gz"),os.path.join(maf_dir,chrN+".maf"))
            rm_cmd="rm {0}".format(os.path.join(maf_dir,chrN+".maf"))
            cmd_list=[gunzip_cmd,] + cmd_list + [rm_cmd]
        combined_cmd_list+=cmd_list
    return combined_cmd_list