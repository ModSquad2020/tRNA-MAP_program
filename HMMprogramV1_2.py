'''
This program is designed to search proteomes using HMMsearch from HMMER
There exists two refined domain profile-HMMs required for this program to work for 
both bacteria and single-cell euekaryotes. 
Please remember to update file_to_open path in def runHmmer!!!
#Still need to review and implement evalue cutoff
<<<<<<< HEAD
Arguements
=======

Arguements

>>>>>>> 04a53dfe53860043c3f3770683bb472b2827ea13
hmmsearch <hmmmodel> <proteome> 
'''

import argparse
import subprocess
import pandas as pd


def parseInput():
    """Parse command line input and output files"""
    
    #Argparse setup
    argParser = argparse.ArgumentParser(description = '')
    argParser.add_argument('-k', '--kingdom', required = True, choices = {'E', 'A', 'B'}, help = 'Domain of organism')
    argParser.add_argument('-p', '--Proteome_input', required = True, help = 'Proteome file input')
    argParser.add_argument('-o', '--output_file', required = True, help = 'Output file name')
    argParser.add_argument('-e', '--Evalue_cutoff', required = False, default = 0.0001, help = 'Define new evalue cutoff ,sci notation can be used') 
    
    #example input
    #python3 HMMprogram.py -k b -p Ecoli_proteome.faa -o ./HMMprofileResults/EcoliHMMsearchResults
    
    clArgs = argParser.parse_args()

    refKingdomType = clArgs.kingdom
    refProteome = clArgs.Proteome_input
    outFile = clArgs.output_file 
    eValue = clArgs.Evalue_cutoff

    return refKingdomType, refProteome, outFile, eValue

def runHMMer(refKingdomType, refProteome):
    """take in correct reference Profile-HMM and Proteome and run HMMsearch"""
    
    file_to_open = None
    if refKingdomType == 'Bacteria':
        file_to_open = './BactProfileHMM.txt' #updatepaths 
    elif refKingdomType == 'Euaryote':
        file_to_open = './EukProfileHMM.txt' #updatepaths
    else:
        print('no set of proteins for archaea yet, sorry!')
    
    try:
        Proteome = refProteome 
        TableOutF = "TableOut" 
        command ="hmmsearch "

        command_string = command + "--tblout " + TableOutF + " " + file_to_open + " " + Proteome

        inresult = subprocess.check_output(command_string, shell=True)

        result = inresult.decode("utf-8")
        return result, TableOutF
        
    except TypeError:
        pass


def readHMMerOutput(TableOutF, eValue):
    '''read table output from HMMer and take domain hits with set evalue and make a dictionary'''
    Hmm_search_dict = {}
    with open(TableOutF, 'r') as TblOut:
        
        table = TblOut.readlines()[3:-10]

        for line in table:
            sLine = line.split()

            if float(sLine[4]) < float(eValue):
                Hmm_search_dict[sLine[2]] = sLine[3:5]
    print(Hmm_search_dict)
    return Hmm_search_dict


    #print(Hmm_search_dict)
    


def domain_finder(Hmm_search_dict, refKingdomType):
    '''make domain dictionary with domain as key'''

    file_to_open = None
    
    if refKingdomType == 'Bacteria':
        file_to_open = './DomainDicBact.txt' #updatepaths 
    elif refKingdomType == 'Eukaryote':
        file_to_open = './DomainDicEuk.txt' #updatepaths
    else:
        print('no set of proteins for archaea yet, sorry!')
        
    try:
        dataFrame = pd.read_csv(file_to_open,sep='\t')

        df = pd.read_csv(file_to_open,sep='\t')
        domain_dict = df.set_index('Domain').T.to_dict('list')
        
        
        return domain_dict

    except TypeError:
        pass


    dataFrame = pd.read_csv(file_to_open,sep='\t')

    df = pd.read_csv(file_to_open,sep='\t')
    domain_dict = df.set_index('Domain').T.to_dict('list')

    #print(domain_dict)
    return domain_dict

'''
    domain_dict = {}
    for key in Hmm_search_dict.keys():
        sub_df = dataFrame[dataFrame["Domain"]==key]

        if sub_df.values.tolist():
            domain_dict[key] = sub_df.values.tolist()


        if sub_df.values.tolist():
            domain_dict[key] = sub_df.values.tolist()

    #print(domain_dict)
    return domain_dict
'''
def combine_dicts(domain_dict, Hmm_search_dict):
    '''read dictionary and merge interesecting domain from 
    tbleout info and domaindictionaries'''
    for domain in domain_dict: 
        if domain in Hmm_search_dict:
            for info in domain_dict[domain]:
                Hmm_search_dict[domain].append(info)           

    #print(Hmm_search_dict)
    return Hmm_search_dict 

def dataframe_out(Hmm_search_dict, outFile):

    outFrame = pd.DataFrame.from_dict(Hmm_search_dict, orient='index')
    rawOutFrame = outFrame.to_csv(outFile, sep='\t',index=True,header=False)
    #print(outFrame)
    return rawOutFrame

'''
def Parse_dataframe(outFile):
    infoDict = {}

    with open('outFileTEST.txt', 'r') as inF:
        lines = inF.readlines()[1:]
        for line in lines:
            line_array = line.split('\t')
            split_info = line_array[2].split(',')
            print(split_info)
            # overwrite the old data
            line_array[2] = split_info[0]
            for i in range(1, len(split_info)):
                line_array.append(split_info[i])

            with open(outFile + ".txt", "a") as output:
                for item in line_array:
                    output.write(item)
                    output.write("\t")
    print(split_info)
  '''              

def main(inCL = True, refKingdomType = None, refProteome = None, outFile = None, eValue = 0.0001):
    
    #Parse command line or recieve args
    
    if inCL == True:
        refKingdomType, refProteome, outFile, eValue = parseInput()

    #run HMMer on sequences
    result, TableOutF = runHMMer(refKingdomType, refProteome)
    
    #Parse HMMer output
    Hmm_search_dict = readHMMerOutput(TableOutF, eValue)
    
    #Assemble dictionary of identified domains
    domain_dict = domain_finder(Hmm_search_dict, refKingdomType)
    
    results_dict = combine_dicts(domain_dict, Hmm_search_dict)
    
    dataframe_out(results_dict, outFile)
    
    return outFile



