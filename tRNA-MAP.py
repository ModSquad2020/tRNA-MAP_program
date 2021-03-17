"""
Call different aspects of tRNA-MAP program


"""

def parseInput():
    """Parse command line input and output files"""
    import argparse
    
    #Argparse setup
    argParser = argparse.ArgumentParser(description = 'Predict tRNA Modifications')
    argParser.add_argument('-s', '--tRNAscan_ss', required = True, help = 'tRNAscan-SE secondary structure output')
    argParser.add_argument('-p', '--protein_sequences', required = True, help = 'protein fasta')
    argParser.add_argument('-o', '--output_directory', required = True, help = 'Path to output results')
    argParser.add_argument('-k', '--kingdom', required = True, choices = {'E', 'A', 'B'}, help = 'Domain of organism')
    argParser.add_argument('--cm_database', required = False, default = './CMlibrary.txt', help = 'File containing information about tRNA mod CMs to search')
    argParser.add_argument('--cp_database', required = False, default = './probLibrary.txt', help = 'File containing information about conditional probabilities')
    argParser.add_argument('--skip_sprinzl_align', required = False, action = 'store_true', help = 'Skips sprinzl alignment step to save time if this has already been done.')
    argParser.add_argument('--e_cutoff', required = False, default = 0.0001, help = 'Cutoff e value for protein searches, default is 0.0001')
    
    kingdoms = {'E': 'Eukaryota', 'A': 'Archaea', 'B': 'Bacteria'}
    #'-s ../SlicerV2/data/secStruct/strePneu_TIGR4-tRNAs.ss.sort -o ./strepneumo_test/ -k B -p ./strepneumo_test/GCF_000006885.1_ASM688v1_protein.faa --skip_sprinzl_align'.split()
    clArgs = argParser.parse_args()
    tRNAstruct = clArgs.tRNAscan_ss
    protSeqs = clArgs.protein_sequences
    orgKing = kingdoms[clArgs.kingdom]
    outDir = clArgs.output_directory
    cmFile = clArgs.cm_database
    probFile = clArgs.cp_database
    skipSprinzl = clArgs.skip_sprinzl_align
    eCutoff = float(clArgs.e_cutoff)
    
    return tRNAstruct, protSeqs, orgKing, outDir, cmFile, probFile, skipSprinzl, eCutoff

def main():
    """"""
    import mapSeqs
    import HMMprogramV1_1
    import predVis
    
    secStruct, protSeqs, orgKing, outDir, cmFile, probFile, skipSprinzl, eCutoff = parseInput()
    
    #format outDir convention
    if outDir[-1] == '/':
        outDir = outDir[:-1]
    
    #Predict tRNA sequences
    seqFile, scoreFile = mapSeqs.main(inCL = False, ssFile = secStruct, domain = orgKing, 
                                      outputFile = outDir, cmsToSearch = cmFile, 
                                      condProbs = probFile, skipSpr = skipSprinzl)
    
    #Search protein sequences
    protFile = '{0}/domainHits.txt'.format(outDir)
    domainHits = HMMprogramV1_2.main(inCL = False, refKingdomType = orgKing, 
                                     refProteome = protSeqs, outFile = protFile, eValue = eCutoff)
    
    domainHits = str()
    predVis.main(inCL = False, calls = seqFile, 
                 scores = scoreFile, proteins = domainHits, 
                 outFile = outDir)
    
    
    
    
main()
