"""
Sort sequences and run CMs; generate output file

1) filter out sequences with wrong base
--> Run sprinzl positions; 
--> Output seqs with correct base to fasta
2) Make predictions using cmsearch
--> Parse cmsearch results from +/- files
3) Sort cmsearch results into output file
--> Write previously-specified output file OR matt's file format

New edit: Parse through search fasta to find tRNA sequences that did not have any hits in the search results. 
 - Make sure to account for them in some way

"""

def parseInput():
    """Parse command line input and output files"""
    import argparse
    
    #Argparse setup
    argParser = argparse.ArgumentParser(description = 'This program generates track hubs that display Modomics data')
    argParser.add_argument('-s', '--tRNAscan_ss', required = True, help = 'list of tRNAscan-SE secondary structure file paths')
    argParser.add_argument('-o', '--output_directory', required = True, help = 'Path to output results')
    argParser.add_argument('-k', '--kingdom', required = True, choices = {'E', 'A', 'B'}, help = 'Domain of organism')
    argParser.add_argument('--cm_database', required = False, default = './CMlibrary.txt', help = 'File containing information about tRNA mod CMs to search')
    argParser.add_argument('--prob_database', required = False, default = './probLibrary.txt', help = 'File containing information about tRNA mod CMs to search')
    argParser.add_argument('--skip_sprinzl_align', required = False, action = 'store_true', help = 'Skips sprinzl alignment step to save time if this has already been done.')
    
    kingdoms = {'E': 'Eukaryota', 'A': 'Archaea', 'B': 'Bacteria'}
    
    #'-s ../SlicerV2/data/secStruct/strePneu_TIGR4-tRNAs.ss.sort -o ./strepneumo_predictions -k B'.split()
    clArgs = argParser.parse_args('-s ../SlicerV2/data/secStruct/strePneu_TIGR4-tRNAs.ss.sort -o ./strepneumo_test -k B'.split())
    tRNAstruct = clArgs.tRNAscan_ss
    orgKing = kingdoms[clArgs.kingdom]
    outDir = clArgs.output_directory
    cmFile = clArgs.cm_database
    probFile = clArgs.prob_database
    skipSprinzl = clArgs.skip_sprinzl_align
    
    return tRNAstruct, orgKing, outDir, cmFile, probFile, skipSprinzl

def sortPositions(inputFile):
    """Read tRNA sprinzl position file"""
    
    sequence = str() #Store tRNA sequence
    posDict = dict()  #Store tRNA sequence as a list
    position = [] #Store sprinzl positions
    
    with open(inputFile, 'r') as inF:
        lines = inF.readlines()
        
        for line in lines[2:]:
            splitLine = line.strip().split('\t')
            
            #Remove spaces for output fasta
            if splitLine[1] != '-':
                sequence += splitLine[1] #Full sequence
            
            posDict[splitLine[0]] = splitLine[1] #Dictionary of positions
            position.append(splitLine[0]) #List of sprinzl positions
        
        inF.close()
    
    return sequence, position, posDict

def writeFasta(outFileName, seqDict, blockSize = 60):
    """Write fasta ouput of unmodified tRNA sequences to use with tRNAscan-SE secondary structure"""
    
    with open(outFileName, 'w') as outF:
        
        for seqName in seqDict.keys():
            
            outF.write('>{0}\n'.format(seqName))
            
            for block in range(0, len(seqDict[seqName]), blockSize):
                try:
                    outF.write('{0}\n'.format(seqDict[seqName][block:block+60]))
                except IndexError:
                    outF.write('{0}\n'.format(seqDict[seqName][block:]))
        
        outF.close()

def cmSearchParser(inputFile):
    """Parse cmsearch --tblout file for e value and sequence name"""
    
    resultsDict = dict()
    
    with open(inputFile, 'r') as inF:
        
        lines = inF.readlines()
        
        #Iterate through lines in file
        for line in lines:
            
            #Skip though metadata
            if line[0] == '#':
                pass
            else:
                splitLine = line.strip().split(' ')
                
                #Remove white space
                values = []
                for x in splitLine:
                    if x != '':
                        values.append(x)
                    else:
                        pass
                
                #Name useful info
                name = values[0]
                cm = values[1]
                score = values[14]
                E = values[15]
                
                resultsDict[name] = score
        
        inF.close()
    
    return resultsDict   

def parseIsoacceptors(inFile):
    """Parse input file type for sorting tRNAs by modification state"""
    """Make sure to remove tRNAs from negative set if it appears in both + and - sets"""
    """Add option to specify specific slices in the tRNA sequence!!!"""
    #Store data
    
    modomicsData = dict()
    
    with open(inFile, 'r') as inF:
        
        entries = inF.read().split('>')[1:] #Split on delimiter and skip empty first string
        
        for entry in entries:
            
            lines = entry.split('\n')[:-1]
            
            #Parse info
            info = lines[0].split('\t')
            
            modomicsData[' '.join(info[0:2])] = dict() #Make a dictionary key for each modified position in file
            
            #Store unmodified base and mod position
            modomicsData[' '.join(info[0:2])]['refBase'] = info[2]
            modomicsData[' '.join(info[0:2])]['refPos'] = info[1]
        
            for line in lines[1:]:
                splitLine = line.split('\t')
                
                #Additional sorting criteria to make parsing easier later on
                species = splitLine[0].split('-')[0]
                isoacceptor = '-'.join(splitLine[0].split('-')[1:4])
                
                #Sort entries in sets
                try:
                    try:
                        modomicsData[' '.join(info[0:2])][species][splitLine[1]].add(isoacceptor)
                    except KeyError:
                        modomicsData[' '.join(info[0:2])][species] = {'+': set(), '-': set(), 'u': set()}
                        modomicsData[' '.join(info[0:2])][species][splitLine[1]].add(isoacceptor)
                except KeyError:
                    modomicsData[' '.join(info[0:2])][species] = {'+': set(), '-': set(), 'u': set()}
                    modomicsData[' '.join(info[0:2])][species][splitLine[1]].add(isoacceptor)
                    
        inF.close()
    
    return modomicsData   

def parseSummary(inputFile):
    """Parse input summary file for CMs"""
    
    infoDict = dict()
    modifiedPositions = dict()
    
    with open(inputFile, 'r') as inF:
        lines = inF.readlines()[1:]
        
        for line, num in zip(lines, range(0, len(lines))):
            
            splitLine = line.strip().split('\t')
            
            infoDict['_'.join([splitLine[0], splitLine[1], splitLine[2]])] = {
                             'mod': splitLine[0], 
                             'pos': splitLine[1], 
                             'refBase': splitLine[2], 
                             'kingdom': splitLine[3], 
                             'posCM': splitLine[4], 
                             'negCM': splitLine[5]}
            try:
                modifiedPositions['_'.join([splitLine[1], splitLine[2]])] += '_'.join([splitLine[0], splitLine[1], splitLine[2]])
            except KeyError:
                modifiedPositions['_'.join([splitLine[1], splitLine[2]])] = ['_'.join([splitLine[0], splitLine[1], splitLine[2]])]
            
        inF.close()
    
    return infoDict, modifiedPositions

def nameIsodecoders(inputFile):
    """
    Read in tRNA isodecoder names and output dictionary of named isodecoders with their sequences and corresponding gene names
    Note: WONT HANDLE INTRONS!!!!!
    """
    
    acceptorScores = {}
    acceptorNames = {}
    
    with open(inputFile, 'r') as inF:
        
        entries = inF.read().split('\n\n')
        
        for entry in entries[:-1]:
            tRNAinfo = entry.split('\n')
            
            #tRNAscan ID name
            tRNAscanID = tRNAinfo[0].split(' ')[0]
            
            #Name isoacceptor
            AA = tRNAinfo[1].split('\t')[0].split(' ')[-1] #Amino acid
            AC = tRNAinfo[1].split(' ')[2] #Anticodon
            isoacceptorName = '-'.join([AA, AC])
            
            #Score
            score = float(tRNAinfo[1].split(' ')[-1])
            
            #sequence
            sequence = tRNAinfo[-2].split(' ')[-1]
            
            try:
                try:
                    acceptorScores[isoacceptorName][sequence] = score
                    acceptorNames[isoacceptorName][sequence].append(tRNAscanID)
                except KeyError:
                    acceptorScores[isoacceptorName][sequence] = score
                    acceptorNames[isoacceptorName][sequence] = [tRNAscanID]
            except KeyError:
                acceptorScores[isoacceptorName] = {sequence: score}
                acceptorNames[isoacceptorName] = {sequence: [tRNAscanID]}
    
    isoacceptors = {}
    isoFastaDict = {}
    scanNameToIsoName = {}
    
    #Sort isodecoders by score and name them
    for acc in acceptorScores.keys():
        for tRNA, pos in zip(acceptorScores[acc].keys(), range(1, len(acceptorScores[acc].keys())+1, 1)):
            
            for score in sorted(acceptorScores[acc].values()):
                
                if acceptorScores[acc][tRNA] == score:
                    tRNAnames = [str(x) for x in acceptorNames[acc][tRNA]]
                    faHeader = '{0}\tScore: {1}\tGene names: {2}'.format('-'.join([str(acc), str(pos)]), str(score), ', '.join(tRNAnames))
                    isoName = '-'.join([acc, str(pos)])
                    isoacceptors[isoName] = {'seq': tRNA, 
                                             'gene names': acceptorNames[acc][tRNA], 
                                             'score': score}
                    isoFastaDict[tRNA.upper()] = isoName
                    
                    for scanName in acceptorNames[acc][tRNA]:
                        scanNameToIsoName[scanName] = isoName
                    
                else:
                    pass
    
    return isoacceptors, isoFastaDict, scanNameToIsoName
    


def parseProbs(inputFile):
    """Parse conditional probability input file"""
    
    #Store conditional probabilities
    probsDict = {}
    
    with open(inputFile, 'r') as inF:
        lines = inF.readlines()
        
        for line in lines[1:]:
            splitLine = line.strip().split('\t')
            
            modLabel = '_'.join(splitLine[0:3])
            residue = splitLine[3]
            pos = splitLine[4]
            posProb = float(splitLine[5])
            negProb = float(splitLine[6])
            
            #Compile dictionary
            try:
                try:
                    probsDict[modLabel][pos][residue] = [posProb, negProb]
                except KeyError:
                    probsDict[modLabel][pos] = {residue: [posProb, negProb]}
                        
            except KeyError:
                probsDict[modLabel] = {pos: {residue: [posProb, negProb]}}
                
        inF.close()
    
    return probsDict
    
def assignProb(tRNAseqs, modInfo, probDict):
    """Calculate conditional probabilities for anticodon modifications"""
    import numpy as np
    mod, refPos, refBase = modInfo.split('_')
    
    probs = {'modScore':{}}
    
    #Iterate through sequences
    for iso in tRNAseqs.keys():
        
        seq = tRNAseqs[iso]
        
        #Iterate through interesting positions
        #for pos in probDict.keys(): ##################################################################
        #Figure out how to combine if multiple positions specified!!!!
        pos = list(probDict.keys())[0]
        
        #Handle U vs T in different datasets
        seqBase = seq[pos]
        if seq[pos] == 'T':
            seqBase = 'U'
        
        try:
            #Get conditional probabilities for observed base
            probs['modScore'][iso] = np.log2(probDict[pos][seqBase][0]/probDict[pos][seqBase][1])
            
        except KeyError:
            if seqBase == '-':
                pass
            else:
                print(pos, seqBase)
          
    return probs
        
def main():
    """Execute commands"""
    import os
    
    #Recieve command line args
    #tRNAstruct, orgKing, outfile, sortFile , runFile
    ssFile, domain, outputFile, cmsToSearch, condProbs, skipSpr = parseInput()
    
    tempDir = '{0}/tmp'.format('/'.join(outputFile.split('/')))
    
    #Directory for output files 
    try:
        os.system('mkdir {0}'.format('/'.join(outputFile.split('/'))))
    except FileExistsError:
        pass
    
    #Make temporary directory for sprinzl position output
    try:
        os.system('mkdir {0}'.format(tempDir))
    except FileExistsError:
        pass
    
    #Parse CM summary file
    allCMs, modInfo = parseSummary(cmsToSearch)
    
    #Parse conditional probability file
    probsDict = parseProbs(condProbs)
    #Add to modInfo
    for modTag in probsDict.keys():
        try:
            modInfo['_'.join(modTag.split('_')[1:])].append(modTag)
        except KeyError:
            modInfo['_'.join(modTag.split('_')[1:])] = [modTag]
    
    #Store tRNA sequences
    tRNAseq = dict() #Store tRNA sequences {tRNAscan-SE ID: sequence}
    tRNApos = dict() #Store tRNA positions {tRNAscan-SE ID:{sprinzl position: base}}
    isoDict = dict() #Store isodecoder sequences and respective tRNAscan-SE IDs {sequence: [tRNAscan-SE ID]}
    
    #Store prediction statistics
    predStats = dict()
    
    #Store tRNA isodecoder names
    faDict = {}
    isoPosDict = {}
    
    isoInfo, fastaSeqs, nameMap = nameIsodecoders(ssFile)
        
    outD = '{0}/sprinzl_alignments'.format(tempDir)
        
    if skipSpr == False:
        try:
            os.system('mkdir {0}'.format(outD))
            os.system('./tRNA_sprinzl_pos -c ./map-sprinzl-pos.conf -d {0} -s {1} -o {2}'.format(domain, ssFile, outD))
        except FileExistsError:
            os.system('./tRNA_sprinzl_pos -c ./map-sprinzl-pos.conf -d {0} -s {1} -o {2}'.format(domain, ssFile, outD))

    #Read through files from sprinzl position namer.
    sprinzlFiles = os.listdir(outD)
        
    #Store CM results for each mod
    predResults = {}
        
    #Iterate through specified CMs:
    for cm in allCMs.keys():
        
        refBase = allCMs[cm]['refBase']
        refPos = allCMs[cm]['pos']
        refMod = allCMs[cm]['mod']
        
        predResults['_'.join([refMod, refPos, refBase])] = {}
        
        #Iterate through sprinzl alignments
        for file in sprinzlFiles:

            #Define relevant file names
            fName = '{0}/{1}'.format(outD, file)
            tRNAscanID = '.'.join(file.split('.')[:-1])

            #Skip tmp files
            if fName.split('.')[-1] == 'pos':

                #Read file
                sequence, position, posDict = sortPositions(fName)

                if posDict[refPos] == refBase:
                    tRNAseq[tRNAscanID] = sequence
                    tRNApos[tRNAscanID] = posDict

                    #Sort by isodecoders:
                    try:
                        isoDict[sequence].append(tRNAscanID)
                    except KeyError:
                        isoDict[sequence] = [tRNAscanID]
                else:
                    pass
            else:
                pass
        
        #Re-pair sequences with isoacceptor names
        for seq in isoDict.keys():
            for ID in isoDict[seq]:
                faDict[nameMap[ID.split('-')[0]]] = seq
                isoPosDict[nameMap[ID.split('-')[0]]] = tRNApos[ID]
        
        #Write modifiable sequences into fasta file for cmSearch
        faName = '{0}/{1}_{2}.fasta'.format(tempDir, refPos, refBase)
        writeFasta(faName, faDict)
            
        #Search sequences against fasta file
        os.system('cmsearch -g --tblout {0}-{1}-hits.txt {2} {0}'.format(faName, allCMs[cm]['posCM'].split('/')[-1], allCMs[cm]['posCM']))
        os.system('cmsearch -g --tblout {0}-{1}-hits.txt {2} {0}'.format(faName, allCMs[cm]['negCM'].split('/')[-1], allCMs[cm]['negCM']))
            
        #Handle when there is no results file (in the case that one CM has no sequences in it)
        try:
            predResults['_'.join([refMod, refPos, refBase])]['posResults'] = cmSearchParser('{0}-{1}-hits.txt'.format(faName, allCMs[cm]['posCM'].split('/')[-1]))
        except FileNotFoundError:
            print('Warning: positive search results file not found')
            posResults = dict()
                
        try:
            predResults['_'.join([refMod, refPos, refBase])]['negResults'] = cmSearchParser('{0}-{1}-hits.txt'.format(faName, allCMs[cm]['negCM'].split('/')[-1]))
        except FileNotFoundError:
            print('Warning: negative search results file not found')
            negResults = dict()
        
        #Add anticodon-loop conditional probabilities
        
        for modTag in probsDict.keys():            
            predResults[modTag] = assignProb(isoPosDict, modTag, probsDict[modTag])
        
    
    #Output file labelling
            
    callFile = '{0}/modCalls.txt'.format(outputFile)
    scoreFile = '{0}/modScores.txt'.format(outputFile)
    
    #write output file
    writeOutput(callFile, scoreFile, predResults, tRNApos, allCMs, modInfo, nameMap)

def writeOutput(callsFile, scoreFile, preds, isos, searchedMods, modPositions, nameMap):
    """Sort isodecoders based on MODOMICS input file"""
    from tRNAinfo import sprinzl2coords as s2c
    
    with open(callsFile, 'w') as cF, open(scoreFile, 'w') as sF:
        
        sprinzlPos = list(s2c.keys())[1:]
        
        cF.write('isoacceptor\t{0}\n'.format('\t'.join(sprinzlPos)))
        sF.write('isoacceptor\t{0}\n'.format('\t'.join(sprinzlPos)))
        
        #reorder into isodecoder names
        isoAlign = {}
        for tRNA in isos.keys():
            isoAlign[nameMap[tRNA.split('-')[0]]] = isos[tRNA]
        
        #Write results to file
        for iso in sorted(isoAlign.keys()):
            cF.write('tRNA-{0}'.format(iso))
            sF.write('tRNA-{0}'.format(iso))
            
            #Warn about abnormal length tRNAs:
            if len(sprinzlPos) != len(isoAlign[iso].keys()):
                print('Warning: tRNA-{0} is not the expected length; skipping gapped regions'.format(iso))
            
            #Iterate through sprinzl positions
            for pos in sprinzlPos:
                
                base = isoAlign[iso][pos]
                
                #calc mod odds score, call mod/unmod based on odds score
                try:
                    
                    modScores = [] #Store multiple scores if present
                    modList = []
                    for mod in modPositions['_'.join([pos, base])]:
                        modScore = float()
                        #calculate mod score
                        if 'modScore' in preds[mod].keys():
                            modScore = preds[mod]['modScore'][iso]
                        else:
                            modScore = float(preds[mod]['posResults'][iso]) - float(preds[mod]['negResults'][iso])
                        
                        #get short name from sorting string
                        shortName = mod.split('_')[0]
                        
                        #Store multiple mod scores
                        modScores.append(modScore)
                        modList.append(shortName)
                   
                    ###################################################
                    #Handle when multiple mods can occur at the same position -> Select the mod that has the highest score
                    for mod, score in zip(modList, modScores):
                        
                        if score == max(modScores):
                            
                            if score > 0:
                                
                                cF.write('\t{0}'.format(mod))
                                sF.write('\t{0}'.format(score))

                            else:
                                cF.write('\t{0}'.format(base))
                                sF.write('\t{0}'.format(score))
                        else:
                            pass 
                    
                    
                    
                    #Write correct mod call
                    #if modScore > 0:
                     #   cF.write('\t{0}'.format(shortName))
                      #  sF.write('\t{0}'.format(modScore))
                            
                    #else:
                     #   cF.write('\t{0}'.format(base))
                      #  sF.write('\t{0}'.format(modScore))
                        
                
                #Write reference base if no modification available
                except KeyError:
                    cF.write('\t{0}'.format(base)) #Write base at unmodified position
                    sF.write('\t-') #Denote no mod prediction at this spot
                    
            cF.write('\n')
            sF.write('\n')
                
        cF.close()
        sF.close()
    
main()
