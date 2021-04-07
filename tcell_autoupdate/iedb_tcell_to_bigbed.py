# This is a script file to download T cell eptiopes from IEDB,
# extract necessary rows and columns and create a beddetail and bigbed
# file for uploading to the UCSC genome browser
# Author: Santrupti Nerli
# email: snerli@ucsc.edu
# Date: March 24, 2021

# load necessary libraries
import os
import re
import sys
import ssl
import wget
import zipfile
import argparse

'''
Command to run this script (I have tested it with python3):

python3 iedb_tcell_to_bigbed.py -path <path where T cell epitopes from IEDB are downloaded and bigbed file is created> -gb_tools <path to gb tools>

'''

# method to download the T cell eptiopes from IEDB webserver
def download_tcell_epitopes(dir):

    # All the t-cell epitopes
    # are available at: https://www.iedb.org/downloader.php?file_name=doc/tcell_full_v3.zip
    # download this file and filter only the t-cell epitopes relevant for SARS-CoV-2
    print ("Downloading T cell epitope file from IEDB with wget module...")

    zip = "tcell_epitopes_iedb.zip"
    target_file = dir+'/'+zip
    url = "https://www.iedb.org/downloader.php?file_name=doc/tcell_full_v3.zip"
    # If we dont use this we will get a SSL certificate verification error
    ssl._create_default_https_context = ssl._create_unverified_context
    # donwload the zip file
    wget.download(url, target_file)

    print ("Finished downloading successfully. Extracting...")

    # extract the zipped file
    with zipfile.ZipFile(target_file,"r") as zipfilehandler:
        zipfilehandler.extractall(dir)

    print ("Finished extracting.")

# method to get input arguments
# typically, the csv is fixed to tcell_full_v3
# as it is in the webserver
# pid is to get the chromosome start and end points
# please use the default file shipped with this script
def get_args():

    """
    Parse command line arguments
    """

    parser = argparse.ArgumentParser(description="Method download T cell epitopes from IEDB and convert to bigbed file")
    parser.add_argument("-path", help="path where csv file will be downloaded")
    parser.add_argument("-csv", help="csv file name of the T cell epitope file donwloaded from IEDB", default="tcell_full_v3.csv")
    parser.add_argument("-hgtable", help="hgtable downloaded from genome browser", default="hgTables.txt")
    parser.add_argument("-gb_tools", help="path to gb_tools", default="./")
    args = parser.parse_args()

    return args

# method to read hgtable obtained from genome browser,
# on the browser, go to Tools,
# select clade: Viruses, genome: SARS-CoV-2, group: UniProt Protein Annotations
# track: Precurs. Proteins, output format: selected fields from primary and related tables
# then click get output
# from a new form, select chrom, chromStart, chromEnd, uniprotName and refSeqProt and
# then click get output
# currently, there are discrepancies between chromosome end numbers between hgTable
# and NCBI entry. So, I download this table and update it so that it matches
# with that of the NCBI entry.
def read_hg_table(args):

    chromStart = {}
    chromEnd = {}
    prot_name = {}
    sequence = {}
    inputfilehandler = open(args.hgtable, 'r')
    for line in inputfilehandler:
        line = line.rstrip()
        fields = line.split()
        if len(fields) > 0:
            if "#" not in line:
                chrom = fields[0]
                refseq = fields[len(fields)-2]
                chromStart[refseq] = fields[1]
                chromEnd[refseq] = fields[2]
                prot_name[refseq] = fields[3:len(fields)-2]
                sequence[refseq] = fields[len(fields)-1]

    inputfilehandler.close()
    return chromStart, chromEnd, prot_name, sequence

# method to read the position of the given peptide
def get_start_pos(peptide, chromStart, prot_seq):

    pid_of_interest = ''
    idx_of_interest = -1
    for pid in prot_seq:
        sequence = prot_seq[pid]
        index = sequence.find(peptide)
        if index != -1:
            idx_of_interest = (index*3) + int(chromStart[pid])
            pid_of_interest = pid
            break

    return pid_of_interest, idx_of_interest

# method to get the last field from a url
# this is required since some of the fields in the
# IEDB column have links with corresponding IDs
# which we want to show on the browser
def get_field_from_url(field):

    fields = field.split('/')
    return fields[len(fields)-1]

# read the downloaded tcell_full_v3.csb file and
# extract what is necessary
def read_csv(csv):

    # dictionary to store the indices of necessary fields
    important_epitope_fields_index = {}

    inputfilehandler = open(csv, 'r')
    header = inputfilehandler.readline()
    fields = header.split(',')
    min_epitope_index = 0
    max_epitope_index = 0
    mhc_name = 0
    assay_type = 0
    # fetch the indices of fields related to epitope, mhc and the assay information
    for i in range(0, len(fields)):
        if fields[i] == "Epitope":
            if min_epitope_index == 0:
                min_epitope_index = i
            max_epitope_index = i
        if fields[i] == "MHC" and mhc_name == 0:
            mhc_name = i
            important_epitope_fields_index['alleleName'] = i
        if fields[i] == "Assay" and assay_type == 0:
            assay_type = i

    # scan for another set of headers to get specific information
    # given as conditions in the if and elif statements
    header = inputfilehandler.readline()
    fields = header.split(',')
    for i in range(min_epitope_index, max_epitope_index+1):
        if "Epitope IRI" in fields[i]:
            important_epitope_fields_index['epitopeID'] = i
        elif "Description" in fields[i]:
            important_epitope_fields_index['epitopeName'] = i
        elif "Object Type" in fields[i]:
            important_epitope_fields_index['epitopeType'] = i
        elif "Starting Position" in fields[i]:
            important_epitope_fields_index['aaStartPos'] = i
        elif "Ending Position" in fields[i]:
            important_epitope_fields_index['aaEndPos'] = i
        elif "Antigen Name" in fields[i]:
            important_epitope_fields_index['antigenName'] = i
        elif "Antigen IRI" in fields[i]:
            important_epitope_fields_index['antigenAccession'] = i
        elif "Parent Protein IRI" in fields[i]:
            important_epitope_fields_index['parentProteinAccession'] = i
        elif "Parent Protein" in fields[i]:
            important_epitope_fields_index['parentProtein'] = i
        elif "Organism Name" in fields[i]:
            important_epitope_fields_index['organismName'] = i
        elif "Parent Species IRI" in fields[i]:
            important_epitope_fields_index['parentOrganismID'] = i
        elif "Parent Species" in fields[i]:
            important_epitope_fields_index['parentOrganism'] = i
        elif "Epitope Comments" in fields[i]:
            important_epitope_fields_index['epitopeComments'] = i

    # get assay type, whether positive, negative or so on for a given row
    i = 0
    while True:
        if "Qualitative Measure" in fields[assay_type+i]:
            important_epitope_fields_index['assayType'] = assay_type+i
            break
        i += 1


    # create a dictionary containing beddetail formatted fields
    # for each epitope from the tcell_full_v3.csv file for
    # SARS-CoV-2 organism
    count = 0
    alleptiopes = {}
    for line in inputfilehandler:
        newline = line.replace('\"', '')
        cleanline = re.sub(r"\s+", ' ', newline)
        fields = cleanline.split(',')
        if "Severe acute respiratory syndrome coronavirus 2" in fields[important_epitope_fields_index['organismName']]:
            attachstring = get_field_from_url(fields[important_epitope_fields_index['epitopeID']]) + '\t' + fields[important_epitope_fields_index['epitopeType']] + '\t' + fields[important_epitope_fields_index['aaStartPos']] + '\t' + fields[important_epitope_fields_index['aaEndPos']] + '\t' + fields[important_epitope_fields_index['antigenName']] + '\t' + get_field_from_url(fields[important_epitope_fields_index['antigenAccession']]) + '\t' + fields[important_epitope_fields_index['parentProtein']] + '\t' + get_field_from_url(fields[important_epitope_fields_index['parentProteinAccession']]) + '\t' + fields[important_epitope_fields_index['organismName']] + '\t' + fields[important_epitope_fields_index['parentOrganism']] + '\t' + get_field_from_url(fields[important_epitope_fields_index['parentOrganismID']]) + '\t' + fields[important_epitope_fields_index['epitopeComments']] + '\t' + fields[important_epitope_fields_index['alleleName']] + '\t' + get_field_from_url(fields[important_epitope_fields_index['assayType']])

            if fields[important_epitope_fields_index['epitopeName']] not in alleptiopes:
                alleptiopes[fields[important_epitope_fields_index['epitopeName']]] = attachstring

    #for a in alleptiopes:
    #    print (a, alleptiopes[a])
    inputfilehandler.close()

    return alleptiopes

# main method
def main(args):

    # download the tcell epitopes from IEDB
    download_tcell_epitopes(args.path)
    # read and filter the tcell epitope csv file
    alleptiopes = read_csv(args.path+'/'+args.csv)
    # get the starting and ending positions of nucleotides, peptides in
    # the SARS-CoV-2 genome
    chromStartList, chromEndList, prot_name, prot_seq = read_hg_table(args)

    # fetch the file tag which will be used
    # to create beddetail, bed and bb files with same name
    outfiletag = args.csv.split('.csv')[0]

    beddetailfilename = outfiletag+'.beddetail'
    bedfilename = outfiletag+'.bed'
    bbfilename = outfiletag+'.bb'

    # remove preexistsing files to avoid confusion
    removefiles = outfiletag+'.b*'
    os.system(f"rm {removefiles}")

    header = []
    beddetailfilehandler = open(beddetailfilename, 'w')

    # write each field in the corresponding beddetail file
    for a in alleptiopes:
        chrom = "NC_045512v2"
        peptide = a
        pid_of_interest, chromStart = get_start_pos(peptide, chromStartList, prot_seq)
        if chromStart != -1:
            chromEnd = str(len(list(peptide))*3+int(chromStart))
            reserved = 0
            name = peptide
            score = 1000
            strand = '+'
            thickStart = str(chromStart)
            thickEnd = str(chromEnd)

            beddetailfilehandler.write(chrom+'\t'+
                    str(chromStart)+'\t'+
                    str(chromEnd)+'\t'+
                    name+'\t'+
                    str(score)+'\t'+
                    strand+'\t'+
                    thickStart+'\t'+
                    thickEnd+'\t'+
                    str(reserved)+'\t'+
                    alleptiopes[a]+"\n")

    beddetailfilehandler.close()

    # use gbtools to convert from beddetail to bed and bigbed
    os.system(f"{args.gb_tools}/bedSort {beddetailfilename} {bedfilename}")
    os.system(f"{args.gb_tools}/bedToBigBed {bedfilename} wuhCor1.sizes {bbfilename} -tab -type=bed9+ -as=iedb.as")

# entry to the program
if __name__ == "__main__":

    args = get_args()
    main(args)
