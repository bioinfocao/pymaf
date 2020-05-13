"""
Convert extracted MAF to fasta format.
Stich MAF blocks into continous fasta sequence.
Each species one line.
"""
# Authors: Huojun Cao <bioinfocao at gmail.com>
# License: BSD 3 clause
import os.path
from sys import exit
import glob
#from bx.align import maf
from mpUtils import getCmdArgs
from mpUtils import mm10_multiz60way_29mammal_species
from mpUtils import hg38_multiz100way_58mammal_species
import traceback
SRC_SPLIT_CHAR = '.'

def printHelp():
    "print help message"
    helpMesage="""
    Convert extracted MAF to fasta format.
    Stich MAF blocks into continous fasta sequence.
    Each species one line. Usage: 
    python mafBlocksToFaSeq.py 
    -i input.maf define a input maf file or directory. For directory this program will loop through all file end with -e (default is .maf)
    -o output_file_dir default is the same as input 
    -s mm10,.../all/multiz60way29mammal/multiz100way58mammal #species to keep in output file, all means keep all species (this is default), multiz60way29mammal means 29mammals, multiz100way58mammal means 58mammals(see below):
    -g yes/no # only keep good MAF files that is at least have two species' seq longer than 60bp 
    -e input file end, default is maf
    -f force fresh if file exist, default is False
    -v 0/1/2 verbose level, default is 1
    Please note that not all species' data will exist for every maf file and output fasta file.
    Example:
    python mafBlocksToFaSeq.py -i NR_038861.1.maf -s all -f True
    """
    print(helpMesage)
    print("multiz60way29mammal: {0}".format(mm10_multiz60way_29mammal_species))
    print("multiz100way58mammal: {0}".format(hg38_multiz100way_58mammal_species))
    exit()

def src_split(src):
    fields = src.split(SRC_SPLIT_CHAR, 1)
    spec = fields.pop(0)
    if fields:
        chrom = fields.pop(0)
    else:
        chrom = spec
    return spec, chrom

def get_species_in_block(block):
    species = []
    for c in block.components:
        spec, chrom = src_split(c.src)
        if spec not in species:
            species.append(spec)
    return species

def get_species_in_maf(input_filename):
    """
    Get all species in maf file
    """
    species = []
    for block in maf.Reader(open(input_filename,'r')):
        for spec in get_species_in_block(block):
            if spec not in species:
                species.append(spec)
    return species

def mafToFasta(input_filename,output_file_dir="",species=[],onlyKeepGoodMaf=True, fresh= False, goodMinSize=60,verbLevel=1):
    """
    Convert maf to fasta
    Each species one line
    """
    species_in_maf_file=get_species_in_maf(input_filename)
    if not species_in_maf_file: return # not good maf file
    if not species:
        species=species_in_maf_file
    species=[s for s in species if s in species_in_maf_file] # keep the order of species the same as MAF file
    if not output_file_dir:
        output_file_dir=os.path.dirname(input_filename)
    headPath, tailFN = os.path.split(input_filename)
    output_filename=os.path.join(output_file_dir, tailFN.replace(".maf",".fasta"))
    if (not fresh) and os.path.isfile(output_filename) and os.path.getsize(output_filename) > 1: return
    # get seq
    if verbLevel>0:
        print("Input maf file: {0}".format(input_filename))
        print("Output fasta file: {0}".format(output_filename))
    if verbLevel>1:
        print("Selected species: {0}".format(species))
    texts = {}
    maf_reader = maf.Reader( open(input_filename,'r') )
    for maf_block in maf_reader:
        for s in species:
            if s not in texts: texts[s]=[]
            c = maf_block.get_component_by_src_start( s ) 
            if c: 
                texts[s].append( c.text )
            else:
                texts[s].append( "-" * maf_block.text_size )
    
    # write to fasta file
    goodCount=0
    if onlyKeepGoodMaf:
        referenceSeq="".join(texts[species[0]])
        minSeqLen=float(len("".join(referenceSeq.split("-"))))*0.1
        for s in species[1:]:
            if texts[s]:
                seq="".join(texts[s])
                seq_rm="".join(seq.split("-"))
                if len(seq)>minSeqLen:
                    goodCount+=1
                if goodCount>=2:
                    break
    if goodCount>=2 or (not onlyKeepGoodMaf):        
        with open(output_filename,'w') as out_file:
            for s in species:
                if texts[s]:
                    out_file.write(">{0}\n{1}\n".format(s,"".join(texts[s])))

if __name__=="__main__":
    
    from sys import argv
    pos_arg_list,named_arg_dict=getCmdArgs(argv)
    if len(pos_arg_list)==1 and len(named_arg_dict)==0:
        printHelp()
    else:
        # get parameters
        if "-i" in named_arg_dict:
            input_file_or_dir=os.path.abspath(named_arg_dict["-i"])
        else:
            printHelp()
        if "-o" in named_arg_dict:
            output_file_dir=named_arg_dict["-o"]
        else:
            if os.path.isdir(input_file_or_dir):
                output_file_dir=input_file_or_dir
            else:
                output_file_dir=os.path.dirname(input_file_or_dir)
        if not os.path.exists(output_file_dir): os.makedirs(output_file_dir)
        if "-g" in named_arg_dict and named_arg_dict["-g"]=="no":
            onlyKeepGoodMaf=False
        else:
            onlyKeepGoodMaf=True
        if "-f" in named_arg_dict and named_arg_dict["-f"]=="True":
            fresh=True
        else:
            fresh=False
        if "-v" in named_arg_dict:
            verbLevel=int(named_arg_dict["-v"])
        else:
            verbLevel=1       
        if "-e" in named_arg_dict:
            fileEnd=named_arg_dict["-e"]
            if fileEnd[0]=="." and len(fileEnd)>1: fileEnd=fileEnd[1:]
        else:
            fileEnd="maf"       
        if "-s" in named_arg_dict:
            if named_arg_dict["-s"]=="multiz60way29mammal":
                species=mm10_multiz60way_29mammal_species
            elif named_arg_dict["-s"]=="multiz100way58mammal":
                species=hg38_multiz100way_58mammal_species
            elif named_arg_dict["-s"]=="all":
                species=[]
            else:
                species=named_arg_dict["-s"].split(",")
        else:
            species=[]
        # execuate mafToFasta
        # if input is directory, loop through all maf file in this directory
        if os.path.isdir(input_file_or_dir):
            input_file_list=glob.glob(os.path.join(input_file_or_dir,"*.{0}".format(fileEnd)))
            for input_file in input_file_list:
                try:
                    mafToFasta(input_file,output_file_dir,species,onlyKeepGoodMaf,verbLevel)
                except Exception as e:
                    traceback.print_exc()      
        # if input is a file        
        else:
            mafToFasta(input_file_or_dir,output_file_dir,species,onlyKeepGoodMaf,fresh,verbLevel)