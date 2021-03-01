'''
This program is designed to search proteomes using HMMsearch from HMMER
There exists two refined domain profile-HMMs required for this program to work for 
both bacteria and single-cell euekaryotes. 
Please remember to update file_to_open path in def runHmmer!!!
#Still need to review and implement evalue cutoff

Arguements

hmmsearch <hmmmodel> <proteome> 
'''

import argparse
import subprocess
import pandas as pd


def parseInput():
    """Parse command line input and output files"""
    
    #Argparse setup
    argParser = argparse.ArgumentParser(description = '')
    argParser.add_argument('-k', '--Kingdom_type', required = True, choices= ("b","e"), help = 'Reference protein domain hmm, b=bacteria, e=eukaryote')
    argParser.add_argument('-p', '--Proteome_input', required = True, help = 'Proteome file input')
    argParser.add_argument('-o', '--output_file', required = True, help = 'Output file name')
    argParser.add_argument('-e', '--Evalue_cutoff', required = False, default = 0.0001, help = 'Define new evalue cutoff ,sci notation can be used') 

    #example input
    #python3 HMMprogram.py -k b -p Ecoli_proteome.faa -o ./HMMprofileResults/EcoliHMMsearchResults

    clArgs = argParser.parse_args()

    refKingdomType = clArgs.Kingdom_type
    refProteome = clArgs.Proteome_input
    outFile = clArgs.output_file
    eValue = clArgs.Evalue_cutoff

    return refKingdomType, refProteome, outFile, eValue

def runHMMer(refProfileHMM, refProteome, outFile):
    """take in correct reference Profile-HMM and Proteome and run HMMsearch"""
    if refProfileHMM == "b":
        file_to_open = './BactProfileHMM.txt' #UPDATE PATH
    else:
        file_to_open = './EukProfileHMM.txt' # UPDATE PATH
    
    Proteome = refProteome 
    OutF = outFile 
    command ="hmmsearch "

    command_string = command + "--tblout " + OutF + " " + file_to_open + " " + Proteome

    print(command_string)

    inresult = subprocess.check_output(command_string, shell=True)

    result = inresult.decode("utf-8")
    return result, OutF




def readHMMerOutput(outF, eValue):
    """Read results from hmmer -tblout output file"""
    
    Hmm_search_dict = {}
    
    with open(outF, 'r') as TblOut:
        
        table = TblOut.readlines()[3:-10]
        #print(table)

        for line in table:
            sLine = line.split()

            if float(sLine[4]) < float(eValue):
                Hmm_search_dict[sLine[2]] = sLine[3:5]
            else:
                Hmm_search_dict[sLine[2]] = sLine[3:5]


            #    if float(sLine[4]) < evalue
               #     Hmm_search_dict[sLine[2]] = sLine[3:5]
               # else:
               #     Hmm_search_dict[sLine[2]] = sLine[3:5]

            # sLine is an array of all the words in a given line
            # sLine = ["NP_011778.1", "-", "Kdo", "PF06293.15", "5.3e-08", "30.9"]
            #if sLine[2] in Hmm_search_dict:
                #if float(sLine[4]) > float(Hmm_search_dict[sLine[2]][4]):
                  #  Hmm_search_dict[sLine[2]] = sLine[3:5]
          #  else:
             #   Hmm_search_dict[sLine[2]] = sLine[3:5]

    #print(Hmm_search_dict)
    return Hmm_search_dict


def domain_finder(Hmm_search_dict, refKingdomType):
    if refKingdomType == "b":
        file_to_open = './DomainDicBact.txt' #updatepaths 
    else:
        file_to_open = './DomainDicEuk.txt' #updatepaths

    dataFrame = pd.read_csv(file_to_open,sep='\t')

    domain_dict = {}

    for key in Hmm_search_dict.keys():
        sub_df = dataFrame[dataFrame["Domain"]==key]

        if sub_df.values.tolist():
            domain_dict[key] = sub_df.values.tolist()

    #print(domain_dict)
    return domain_dict

def combine_dicts(domain_dict, Hmm_search_dict):

    missing_domains = []
    for domain in domain_dict: 
        if domain not in Hmm_search_dict:
            missing_domains.append(domain)
        else:
            Hmm_search_dict[domain].append(domain_dict[domain])

    #print(Hmm_search_dict)
    print(missing_domains)
    return missing_domains, Hmm_search_dict 

def dataframe_out(Hmm_search_dict):

    outFrame = pd.DataFrame.from_dict(Hmm_search_dict, orient='index')

    return outFrame.to_csv('output.tsv',sep='\t',index=False,header=False)





if __name__ == "__main__":
    refKingdomType, refProteome, outFile, eValue = parseInput()
    result, OutF = runHMMer(refKingdomType, refProteome, outFile)
    Hmm_search_dict = readHMMerOutput(OutF, eValue)
    domain_dict = domain_finder(Hmm_search_dict, refKingdomType)
    missing_domains, Hmm_search_dict = combine_dicts(domain_dict, Hmm_search_dict)
    dataframe_out(Hmm_search_dict)

    




    




