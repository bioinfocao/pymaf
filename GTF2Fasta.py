"""
Take known or novel transcripts GTF as input,
Return transcripts' fasta sequence extracted from multiple species whole genome alignments MAF file (such as UCSC multiz100way)


Requirements: ucsc-mafsinregion. Download and install as:
http://hgdownload.cse.ucsc.edu/admin/exe/
https://anaconda.org/bioconda/ucsc-mafsinregion
"""
import os
import os.path
from sys import exit
from time import gmtime, strftime
datetime_text = strftime("%Y_%m_%d_%H_%M_%S", gmtime())
import platform
import shutil

from mpUtils import excu, getCmdArgs, getArgsVal, gtfToBed6, catExonsMAFCmds
from mpUtils import exonBedSplitAndChangeName, mafExtractCmds
from mpUtils import mm10_multiz60way_29mammal_species
from mpUtils import hg38_multiz100way_58mammal_species

def printHelp():
    "print help message"
    helpMesage="""Take known or novel transcript GTF as input,
    return transcripts' fasta sequence extracted from whole genome multiple alignments MAF file.
    Usage: python GTF2Fasta.py
    -i inout known/novel GTF file
    -m maf files directory # multiz MAF files dir (chr1.maf.gz ... chrY.maf.gz)
    -n yes/no filterChrNNN # limit to chrN,chrNN,chrNNN. Any chromosome name longer than 6 letters are skipped. default is yes.
    -z yes/no # the MAF file in zipped format (as *.maf.gz this is default) file or in unzipped foramt .maf (no)
    -o output_fasta_dir # the dir for transcripts fasta output, if omitted create extracted_fasta_files folder under input GTF file.
    -t transcript_catMAF_dir # if not specified will put in output_fasta_dir folder
    -f yes/no # make fresh transcript_catMAF file? if no, will skip maf file that already exist (size>1). default is no.
    -s mm10,.../all/multiz60way29mammal/multiz100way58mammal # selected species names, default is keep all species.
       multiz60way29mammal means 29 mammals, multiz100way58mammal means 58 mammals(see below), all means keep all species, which is default.
       Please note that not all species' data will exist for every maf file
    -c transcripts_gtf_to_fasta_cmds.sh # output .sh cmd file, without this argument cmds will be excuate automatically.
    -g yes/no # only keep MAF files that is at least have another two species' fasta sequence longer than 1/10 of the frist(reference) species, default is yes.
    -d yes/no # delete temporary nonoverlapExons_bedFiles_ files, default is yes
    -v 0/1/2 # verbose level default is 1, 0 means nothing, 2 is detail.
    Example:
    python GTF2Fasta.py -i hg38_RefSeq.gtf -m ~/data/refgenomes/multiz100way_hg38 -o ~/data/hg38_transcript_catMAF_fasta/ -c hg38_RefSeq_gtf_to_fasta_cmds.sh
    """
    print(helpMesage)
    print("multiz60way29mammal: {0}".format(mm10_multiz60way_29mammal_species))
    print("multiz100way58mammal: {0}".format(hg38_multiz100way_58mammal_species))
    exit()

def set_mafsInRegion_bin():
    """
    get mafsInRegion bin file
    """
    mafsInRegion_bin=excu("which mafsInRegion", retErr=False, verbose=False)
    if not mafsInRegion_bin:
        if platform.system().find("Darwin")>=0:
            mafsInRegion_bin=os.path.join(code_dir,"bins/mafsInRegion_mac64")
        else:
            mafsInRegion_bin=os.path.join(code_dir,"bins/mafsInRegion_linux64")
    mafsInRegion_bin=mafsInRegion_bin.rstrip()
    return mafsInRegion_bin

def mainFunc(input_file, maf_dir, filterChrNNN, extracted_chr_start_end_maf_dir, output_dir, \
             transcript_catMAF_dir, freshTransCatMAF, cmd_file, species, rmTempBed, gzInput, verbLevel):
    """ The main function
    """
    # find mafsInRegion
    mafsInRegion_bin=set_mafsInRegion_bin()
    # convert gtf to bed6 is input is gtf
    bed_input_file=input_file.replace(".gtf","")+".bed"
    gtfToBed6(input_file,bed_input_file,verbLevel=verbLevel)
    # prepare_nonoverlapExons_BedFiles
    exonBedSplitAndChangeName(bed_input_file, extracted_chr_start_end_maf_dir,temp_nonoverlapExons_bedFiles_dir)
    # make maf extractionn cmds
    maf_extraction_cmd_list=mafExtractCmds(mafsInRegion_bin,temp_nonoverlapExons_bedFiles_dir,\
                                           maf_dir,extracted_chr_start_end_maf_dir,gzInput,verbLevel)
    # make concatenate exons' MAF into transcripts' MAF
    maf_cat_cmd_list=catExonsMAFCmds(bed_input_file,extracted_chr_start_end_maf_dir,transcript_catMAF_dir,freshTransCatMAF,filterChrNNN)
    # convert transcripts' MAF files to fasta files
    maf_to_fasta_cmd="python {0} -i {1} -o {2} -s {3} -g {4} -e maf -v {5}".format(\
        os.path.join(code_dir,"mafBlocksToFaSeq.py"), transcript_catMAF_dir, output_dir, species, onlyKeepGoodMaf, verbLevel)
    # remove temporary files
    rm_tempFolder_cmds=[]
    if rmTempBed: rm_tempFolder_cmds.append("rm -r {0}".format(temp_nonoverlapExons_bedFiles_dir))
    # generate cmd list
    cmd_list=maf_extraction_cmd_list + maf_cat_cmd_list + [maf_to_fasta_cmd,]+ rm_tempFolder_cmds
    if not cmd_file:
        for cmd in cmd_list:
            out,err=excu(cmdText, retErr=True, verbose=True)
            print("Outputs:"+out)
            print("Errors"+err)
    else:
        header_line="#!/bin/sh\n"
        with open(cmd_file,'w') as f:
            f.write(header_line)
            f.write("\n".join(cmd_list))
        print("Please execuate following script file:\n{0}".format(cmd_file))

if __name__=="__main__":
    
    print("Take known or novel transcript GTF as input and return transcripts fasta sequence extracted from MAF file")
    # setting 
    code_dir=os.path.dirname(os.path.abspath(__file__)) #get current script file path
    temp_dir=os.path.join(code_dir, "tempfiles") # temp dir
    if not os.path.exists(temp_dir): os.makedirs(temp_dir)    
    # Get arguments
    from sys import argv
    pos_arg_list,named_arg_dict=getCmdArgs(argv)
    if len(pos_arg_list)==1 and len(named_arg_dict)==0:
        printHelp()
    else:  
        # get parameters
        name_list=["-i", "-m", "-n", "-z", "-o", "-t", "-f", "-s", "-c", "-g", "-d", "-v"]
        input_file, maf_dir, filterChrNNN, gzInput, output_dir, transcript_catMAF_dir, freshTransCatMAF,\
            species, cmd_file, onlyKeepGoodMaf, rmTempBed, verbLevel = getArgsVal(named_arg_dict, name_list)
        for el in [input_file, maf_dir]:
            if el is None:
                printHelp()
        if filterChrNNN and filterChrNNN =="no":
            filterChrNNN=False
        else:
            filterChrNNN=True
        if gzInput and gzInput=="no":
            gzInput=False
        else:
            gzInput=True
        if output_dir is None:
            output_dir=os.path.join(os.path.dirname(input_file),"extracted_fasta_files")
        if not os.path.exists(output_dir): os.makedirs(output_dir)
        if transcript_catMAF_dir is None:
            transcript_catMAF_dir=output_dir
        if not os.path.exists(transcript_catMAF_dir): os.makedirs(transcript_catMAF_dir)
        if freshTransCatMAF and freshTransCatMAF=="yes":
            freshTransCatMAF=True
        else:
            freshTransCatMAF=False
        if species is None: species="all"
        if onlyKeepGoodMaf is None: onlyKeepGoodMaf="yes"
        if rmTempBed and rmTempBed=="no":
            rmTempBed=False
        else:
            rmTempBed=True
        if verbLevel:
            verbLevel=int(verbLevel)
        else:
            verbLevel=1
    ################################################################
    # temp dir for nonoverlapExons_bedFiles
    temp_nonoverlapExons_bedFiles_dir=os.path.join(temp_dir, "nonoverlapExons_bedFiles_{0}".format(datetime_text))
    if not os.path.exists(temp_nonoverlapExons_bedFiles_dir): os.makedirs(temp_nonoverlapExons_bedFiles_dir)
    # extracted_chr_start_end_maf_dir
    extracted_chr_start_end_maf_dir=os.path.join(maf_dir,"extracted_chr_start_end_maf_files")
    if not os.path.exists(extracted_chr_start_end_maf_dir): os.makedirs(extracted_chr_start_end_maf_dir)
    ################################################################
    print("Input GTF file: {0}".format(input_file))
    print("Extract MAF Fasta Seqs to: {0}".format(output_dir))
    mainFunc(input_file, maf_dir, filterChrNNN, extracted_chr_start_end_maf_dir, output_dir, transcript_catMAF_dir, freshTransCatMAF, cmd_file, species, rmTempBed, gzInput, verbLevel)
