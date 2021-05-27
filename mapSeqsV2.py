"""
Improvements:

Must add in inserted regions. to do this:
1) make a dictionary with the regular sprinzl positions that gaps occur aferwatds
2) iterate through this dict and generate a complete list of all the positions in the alignment
3) instead of iterating through the reference list of sprinzl positions in the tRNAinfo.py, iterate through this list while writing output

In Progress:
Adding in protein HMM stuff
- need to make sure code accounts for when two hmms hit to the same protein -> call enzyme as whichever hmm has a better bit score
- For some reason CMs not searching for s4U because fasta is empty.
Linear regression:
- Using onehot encoding because numberical encoding apparently can cause some false correlations
- ['A', 'U'] -> []
- FIX LINEAR REGRESSION SO NUMBER OUTPUTS MAKE SENSE!!!
Output:
- finish writing mod call function. 
 - Change so that all possible predicted mods are written as comma-separated lists
 - Kicks can down the road of predicting which tRNAs have a combined modification (like ms2i6A/ms2t6A)
 - 

"""

def parseInput():
    """Parse command line input and output files"""
    import argparse
    
    #Argparse setup
    argParser = argparse.ArgumentParser(description = 'This program generates track hubs that display Modomics data')
    argParser.add_argument('-s', '--tRNAscan_ss', required = True, help = 'list of tRNAscan-SE secondary structure file paths')
    argParser.add_argument('-o', '--output_directory', required = True, help = 'Path to output results')
    argParser.add_argument('-d', '--domain', required = True, choices = {'E', 'A', 'B'}, help = 'Domain of organism')
    argParser.add_argument('-p', '--proteins', required = True, help = 'Query organism protein sequences')
    argParser.add_argument('--mod_database', required = False, default = './mod_library.txt', help = 'File containing information about tRNA mod CMs to search')
    argParser.add_argument('--skip_sprinzl_align', required = False, action = 'store_true', help = 'Skips sprinzl alignment step to save time if this has already been done.')
    argParser.add_argument('--e_cutoff', required = False, default = '1e-6', help = 'Only accept proteins with an E value cutoff above specified number')
    argParser.add_argument('--cpu', required = False, default = '4', help = 'Number of CPUs to use. Default is 4')
    
    #Store domain names
    domains = {'E': 'Eukaryota', 'A': 'Archaea', 'B': 'Bacteria'}
    
    #'-s ./secStruct/strePneu_TIGR4-tRNAs.ss.sort -o ./validation_testing/strePneu-lmLibrary -p ./test_proteomes/strePneu-GCF_000006885.1_ASM688v1_protein.faa -d B --mod_database ./lm_mod_library.txt'.split()
    #'-s ../SlicerV2/data/secStruct/strePneu_TIGR4-tRNAs.ss.sort -o ./work/strepneumo_test/ -p ./work/strepneumo_test/GCF_000006885.1_ASM688v1_protein.faa -d B --cpu 3'.split()
    #'-s ./secStruct/strePneu_TIGR4-tRNAs.ss.sort -o ./work/strepneumo_test/ -p ./test_proteomes/strePneu-GCF_000006885.1_ASM688v1_protein.faa -d B'.split()
    clArgs = argParser.parse_args()
    tRNAstruct = clArgs.tRNAscan_ss
    orgDom = domains[clArgs.domain.upper()]
    outDir = clArgs.output_directory
    modFile = clArgs.mod_database
    prots = clArgs.proteins
    skipSprinzl = clArgs.skip_sprinzl_align
    eValCutoff = float(clArgs.e_cutoff)
    cpus = clArgs.cpu
    
    return tRNAstruct, orgDom, outDir, modFile, skipSprinzl, eValCutoff, prots, cpus

def bit2prob(bitScore):
    """Convert bit score to probability, assuming a binary outcome"""
    
    return 1/(1+2**(-bitScore))

def lmParse(lmFile):
    """Parse tab-separated linear model file"""
    with open(lmFile, 'r') as inF:
        
        #parsed lm file
        lines =  inF.readlines()
        
        modState = []
        data = []
        
        for line in lines[1:]:
            splitLine = line.strip().split('\t')
            
            modState.append([int(splitLine[1])])
            data.append(splitLine[2:])
        
        #lm information
        mod = lines[0].strip().split('\t')[0]
        columns = lines[0].strip().split('\t')[2:]
        
        return mod, columns, modState, data  

def lmPredict(alignDict, lmFile):
    """Builds linear model from linear-regression training data; searches tRNAs against linear model"""
    import numpy as np
    from sklearn.preprocessing import OneHotEncoder as OHE
    from sklearn.linear_model import LinearRegression as LR
    
    mod, cols, modState, data = lmParse(lmFile)
    
    #Add prediction set to predictor set for encoding.
    #Iterate through and predict tRNAs
    
    trainingCutoff = len(data)
    
    predictionVars = []
    for tRNA in sorted(alignDict.keys()): ################This should be temporary I think
        
        #Collect predictor variables into a list
        tRNAvars = []
        for alignPos in cols:
            try:
                tRNAvars.append(alignDict[tRNA][alignPos])
            except KeyError:
                tRNAvars.append(tRNA.split('-')[0])
        
        predictionVars.append(tRNAvars)
    
    predVars = np.array(data + predictionVars)
    
    predState = np.array(modState)
    
    #Encode each set of elements in the linear regression dataset into onehot 
    encoding = OHE(sparse=False)
    encodedArray = encoding.fit_transform(predVars)
    
    #Instantiate linear model  ###########################NOTE: I MAY HAVE THIS BACKWARDS. NOT EXACTLY SURE WHICH WAY DATA GETS FED INTO LM!!!
    linearModel = LR().fit(encodedArray[:trainingCutoff], predState)
    
    #Linear model statistics
    #print(linearModel.score(encodedArray, predState))
    #print(linearModel.coef_)
    #print(linearModel.intercept_)
    
    predDict = {} #Store tRNA modification lm predictions
    
    for tRNA, pos in zip(sorted(alignDict.keys()), range(0, len(sorted(alignDict.keys())), 1)): #Sorted in the same manner that entries are placed into the list.
        
        #Get specific tRNA data
        tRNAencode = [encodedArray[trainingCutoff + pos]]
          
        #Predict using linear model
        tRNAprob = linearModel.predict(tRNAencode)
        
        predDict[tRNA] = 0.5*(float(tRNAprob)+1) #######THIS IS A PLACEHOLDER!! GOTTA FIGURE OUT HOW TO CONVERT TO PROBABILITIES!!!!!
    
    return predDict

def writeTSV(fName, data, columns, delim = '\t', nd = '-'): #add in columns.
    """Write a delimited text file with a given list of columns with entries as rows; and sub-entries as each data point"""
    
    with open(fName, 'w') as outF:
        
        outF.write(delim.join(columns) + '\n')
        
        for line in sorted(list(data.keys())):
            
            outF.write(line)
            
            for pos in columns[1:]:
                
                #Write data or no data char if col not present
                try:
                    try:
                        
                        outF.write(delim + data[line][pos])
                    except TypeError:
                        try: #Handle when we have a float
                            outF.write(delim + str(round(data[line][pos], 2)))
                        except TypeError: #Handle when you have a protein sequence
                            outF.write(delim + str(round(data[line][pos]['score'], 2)))
                except KeyError:
                    outF.write(delim + nd)   
                
            outF.write('\n')
        
        outF.close()
        
def getPositions(posDict):
    """Build a list of sprinzl positions that includes gaps specific to the organism of interest"""
    from tRNAinfo import sprinzl2coords
    
    #Store gap positions under correct sprinzl pos
    gapsDict = {x: dict() for x in sprinzl2coords.keys()}
    
    #Iterate throught tRNAs to get gap positions
    for tRNA in posDict.keys():
        for pos in posDict[tRNA]:
            
            gapPos = pos.split(':i')[0]
            gapNo = pos.split(':i')[-1]
            #Add gap to dict only if it isn't a sprinzl position
            
            if pos not in gapsDict.keys():
                gapsDict[gapPos][gapNo] = pos
            else:
                pass
    
    #Reformat into list
    allPos = []
    
    for pos in gapsDict.keys():
        allPos.append(pos)
        
        for gap in gapsDict[pos].keys():
            allPos.append(gapsDict[pos][gap])
        
    return allPos

def predictProteins(protSeqs, HMM, eCutoff, temp, cpus):
    """
    1) Use hmmscan to search for tRNA modification enzymes
    2) Parse hmmscan output file
    3) Select best sequence (maybe pick the one with the best % identity or query coverage given the e value is above the cutoff)
    """
    import os
    predFile = '{0}/{1}-hmmResults.txt'.format(temp, '.'.join(HMM.split('/')[-1].split('.')[:-1]))
    
    os.system('hmmsearch --cpu {3} --tblout {0} {2} {1}'.format(predFile, protSeqs, HMM, cpus))
    
    protHits = cmSearchParser(predFile, prot = True)
    
    #Return the best hit or return false if no hit above e value threshold
    hitScores = set( protHits[x][0] for x in protHits.keys())
    for hit in protHits.keys():
        
        if protHits[hit][1] <= eCutoff and protHits[hit][0] == max(hitScores):
            return {'hit': hit, 'score': protHits[hit][2]}
        elif protHits[hit][1] > eCutoff and protHits[hit][0] < max(hitScores):
            return None

        
def searchCMs(allCMs, nameMap, outD, tempDir, cpus):
    """Use cmsearch to generate modification scores. TRANSFER CODE FROM MAIN FUNCTION INTO THIS FUNCTION"""
    
    import os
    
    #Store tRNA sequences
    tRNAseq = dict() #Store tRNA sequences {tRNAscan-SE ID: sequence}
    tRNApos = dict() #Store tRNA positions {tRNAscan-SE ID:{sprinzl position: base}}
    isoDict = dict() #Store isodecoder sequences and respective tRNAscan-SE IDs {sequence: [tRNAscan-SE ID]}
    
    #Store prediction statistics
    predStats = dict()
    
    #Store tRNA isodecoder names
    faDict = {}
    isoPosDict = {}
    
    #Read through files from sprinzl position namer.
    sprinzlFiles = os.listdir(outD)
    
    #Store CM results for each mod
    predResults = {}
    
    #Store positions for each tRNA
    posLists = {}
    
    #Iterate through sprinzl alignments
    for file in sprinzlFiles:

        #Define relevant file names
        fName = '{0}/{1}'.format(outD, file)
        tRNAscanID = '.'.join(file.split('.')[:-1])

        #Skip tmp files
        if fName.split('.')[-1] == 'pos':

            #Read file
            sequence, position, posDict = sortPositions(fName)
            
            tRNAseq[tRNAscanID] = sequence
            tRNApos[tRNAscanID] = posDict
            posLists[tRNAscanID] = position
                
            #Sort by isodecoders:
            try:
                isoDict[sequence].append(tRNAscanID)
            except:
                isoDict[sequence] = [tRNAscanID]
            
        else:
            pass
        
    #Re-pair sequences with isoacceptor names
    for seq in isoDict.keys():
        for ID in isoDict[seq]:
            faDict[nameMap[ID.split('-')[0]]] = seq
            isoPosDict[nameMap[ID.split('-')[0]]] = tRNApos[ID]
    
    #Iterate through specified CMs:
    for cm in allCMs.keys():
        
        refBase = allCMs[cm]['refBase']
        refPos = allCMs[cm]['pos']
        refMod = allCMs[cm]['mod']
        
        predResults['_'.join([refMod, refPos, refBase])] = {}
        
        #Write modifiable sequences into fasta file for cmSearch
        faName = '{0}/{1}_{2}.fasta'.format(tempDir, refPos, refBase)
        writeFasta(faName, faDict)
            
        #Search sequences against fasta file
        os.system('cmsearch -g --cpu {3} --tblout {0}-{1}-hits.txt {2} {0}'.format(faName, allCMs[cm]['posCM'].split('/')[-1], allCMs[cm]['posCM'], cpus))
        os.system('cmsearch -g --cpu {3} --tblout {0}-{1}-hits.txt {2} {0}'.format(faName, allCMs[cm]['negCM'].split('/')[-1], allCMs[cm]['negCM'], cpus))
            
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
        
    return predResults, isoPosDict, tRNApos, posLists

def parseSummary(inputFile):
    """Parse input summary file for CMs"""
    
    cmDict = dict()
    probsDict = dict()
    consDict = dict() # Store constitutive mods
    protDict = dict() # Store the associated protein HMMs for a given mod
    modifiedPositions = dict()
    
    with open(inputFile, 'r') as inF:
        lines = inF.readlines()
        
        for line in lines:
            
            splitLine = line.strip().split('\t')
            
            #Add identifiable labels for the data in each mod
            modTag = '_'.join([splitLine[1], splitLine[2], splitLine[3]])
            posTag = '_'.join([splitLine[2], splitLine[3]])
            
            if splitLine[0].lower() == 'cm': #Parse CM mods
            
                cmDict['_'.join([splitLine[1], splitLine[2], splitLine[3]])] = {
                                 'mod': splitLine[1], 
                                 'pos': splitLine[2], 
                                 'refBase': splitLine[3], 
                                 'posCM': splitLine[4], 
                                 'negCM': splitLine[5]}
                try:
                    modifiedPositions[posTag].add(modTag)
                except KeyError:
                    modifiedPositions[posTag] = set()
                    modifiedPositions[posTag].add(modTag)
                    
                protDict[modTag] = splitLine[6]
                    
                    
            elif splitLine[0].lower() == 'lm': #parse conditional probabilities (change to lin reg soon)
                
                lm = splitLine[4]
                prot = splitLine[5]
                
                #Add protein HMM file
                protDict[modTag] = prot
                
                #Append linear model file
                probsDict[modTag] = lm
                
                try:
                    modifiedPositions[posTag].add(modTag)
                except KeyError:
                    modifiedPositions[posTag] = set()
                    modifiedPositions[posTag].add(modTag)
                    
            elif splitLine[0].lower() == 'cs':
                
                mod = splitLine[1]
                pos = splitLine[2]
                refBase = splitLine[3]
                hmmFile = splitLine[4]
                protDict[modTag] = hmmFile
                
                consDict['_'.join([mod, pos, refBase])] = hmmFile
                
                try:
                    modifiedPositions[posTag].add(modTag)
                except KeyError:
                    modifiedPositions[posTag] = set()
                    modifiedPositions[posTag].add(modTag)
                
            else:
                pass
                
    
        inF.close()
    
    return cmDict, probsDict, consDict, protDict, modifiedPositions
        
def sortPositions(inputFile):
    """Read tRNA sprinzl position file"""
    
    sequence = str() #Store tRNA sequence
    posDict = dict()  #Store tRNA sequence as a list
    position = [] #Store sprinzl positions
    
    with open(inputFile, 'r') as inF:
        lines = inF.readlines()
        
        for line in lines[2:]:
            splitLine = line.strip().split('\t')
            
            #Convert to RNA
            if splitLine[1] == 'T':
                splitLine[1] = 'U'
            
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

def cmSearchParser(inputFile, prot = False):
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
                
                if prot == False:
                    #Name useful info
                    name = values[0]
                    cm = values[1]
                    score = values[14]
                    E = values[15]

                    resultsDict[name] = score
                else:
                    #Name useful info
                    name = values[0]
                    hmm = values[2]
                    score = float(values[5])
                    E = float(values[4])

                    resultsDict[name] = [score, E, score]
        
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
    
def modifySeqs(preds, isos, searhchedMods, modPositions, nameMap, allPos, protHits):
    """Compile mod predictions from different sources"""
    # 1) iterate through the tRNAs in the isodecoder set
    # 2) at each position check if a mod is possible
    #  a) if 2+ mods possible, select one with best score. (how to decide non-probabilistic and probabilistic mods occur together?)
    #  b) if 1 mod possible, see if score is more likely or less likely to happen. Write more likely mod
    # 3) if mod possible, 
    
    tRNAseqs = dict() #Store modified tRNA sequences
    posScores = dict() #Store associated mod scores
    protPos = dict() #Store predicted protein accessions
    protScores = dict() #Store protein hit bit scores
    
    #reorder into isodecoder names
    isoAlign = {}
    for tRNA in isos.keys():
        isoAlign[nameMap[tRNA.split('-')[0]]] = isos[tRNA]
            
    #Write results to file
    for iso in sorted(isoAlign.keys()):
        tRNAseqs[iso] = {}
        posScores[iso] = {}
        protPos[iso] = {}
        protScores[iso] = {}
        
        #Iterate through sprinzl positions in the tRNA
        for pos in allPos:
            
            #Handle extra gaps
            try:
                base = isoAlign[iso][pos].upper()
                
                if base.upper() == 'T':
                    base = 'U'
            
            except KeyError:
                base = '-'
            
            posTag = '_'.join([pos, base])
            
            if posTag in modPositions.keys():
                
                posHits = []
                
                for modTag in modPositions[posTag]:
                    
                    if protHits[modTag] != None:
                        
                        try: #Handle when covariance model is used
                            cmResults = preds['cm']
                            
                            modScore = float()
                            
                            try: #Handle when hits to both are present
                                modScore = float(preds['cm'][modTag]['posResults'][iso])-float(preds['cm'][modTag]['negResults'][iso])
                                
                            except KeyError:
                                
                                try: #No positive hit
                                    modScore = float( 0 - float(preds['cm'][modTag]['negResults'][iso]))
                                    
                                except KeyError:
                                    
                                    try: #no negative hit
                                        modScore = float(float(preds['cm'][modTag]['posResults'][iso]) - 0)
                                        
                                    except KeyError: #No hit to either
                                        
                                        exceptionBreaker = preds['cm'][modTag]
                                        modScore = 0 #############################################Change this do a no data value
                                       
                            posHits.append([modTag, round(bit2prob(modScore), 2), protHits[modTag]['hit'], protHits[modTag]['score']])
                        
                        except KeyError:
                            
                            try: #Handle conditional prob (and soon linear regression)
                                posHits.append([modTag, round(preds['lm'][modTag][iso], 2), protHits[modTag]['hit'], protHits[modTag]['score']])
                                
                            except KeyError: #Handle constitutive modifications
                                posHits.append([modTag, 1, protHits[modTag]['hit'], protHits[modTag]['score']])
                    
                    
                        #Sorted score values
                        sortedHits = sorted(posHits, key = lambda score: (score[1], score[0], score[2], score[3]), reverse = True) 

                        sortedLabels = [x[0].split('_')[0] for x in sortedHits]
                        sortedScores = [str(x[1]) for x in sortedHits]
                        sortedProtes = [x[2] for x in sortedHits]
                        protScoresList = [str(x[3]) for x in sortedHits]

                        #Add prediction data to sequence
                        tRNAseqs[iso][pos] = ','.join(sortedLabels)
                        posScores[iso][pos] = ','.join(sortedScores)
                        protPos[iso][pos] = ','.join(sortedProtes)
                        protScores[iso][pos] = ','.join(protScoresList)
                    
                    
                    else:
                        
                        tRNAseqs[iso][pos] = base
                        posScores[iso][pos] = '-'
                        protPos[iso][pos] = '-'
                        protScores[iso][pos] = '-'
                
            else:
                tRNAseqs[iso][pos] = base
                posScores[iso][pos] = '-'
                protPos[iso][pos] = '-'
                protScores[iso][pos] = '-'
    
    return tRNAseqs, posScores, protPos, protScores
    
def main(inCl = True, ssFile = None, domain = None, 
         outputFile = None, modInfo = './CMlibrary.txt', skipSpr = False, 
         eCut = 1e-6, protSeqs = None):
    """Execute commands"""
    
    import os
    
    #Recieve command line args
    #tRNAstruct, orgKing, outfile, sortFile , runFile
    if inCl == True:
        ssFile, domain, outputFile, cmsToSearch, skipSpr, eCut, protSeqs, cpus = parseInput()
    
    #Remove end slash if present
    if outputFile[-1] == '/':
        outputFile = outputFile[:-1]
    
    tempDir = '{0}/tmp'.format(outputFile)
    
    #Directory for output files 
    try:
        os.system('mkdir {0}'.format(outputFile))
    except FileExistsError:
        pass
    
    #Make temporary directory for sprinzl position output
    try:
        os.system('mkdir {0}'.format(tempDir))
    except FileExistsError:
        pass
    
    #Parse CM summary file
    allCMs, probsDict, constDict, protDict, modInfo = parseSummary(cmsToSearch)
    
    #Give isodecoders GtRNAdb names
    isoInfo, fastaSeqs, nameMap = nameIsodecoders(ssFile)
    
    outD = '{0}/sprinzl_alignments'.format(tempDir)
        
    if skipSpr == False:
        try:
            os.system('mkdir {0}'.format(outD))
            os.system('./tRNA_sprinzl_pos -c ./map-sprinzl-pos.conf -d {0} -s {1} -o {2}'.format(domain, ssFile, outD))
        except FileExistsError:
            os.system('./tRNA_sprinzl_pos -c ./map-sprinzl-pos.conf -d {0} -s {1} -o {2}'.format(domain, ssFile, outD))

    #Run cmsearch on organism's tRNA sequences
    predResults, isoPosDict, tRNApos, positionLists = searchCMs(allCMs, nameMap, outD, tempDir, cpus)
    allPreds = {'cm': predResults, 'lm': {}, 'cs': {}}
    
    #Add anticodon-loop conditional probabilities to the cm prediction dictionary
    for modTag in probsDict.keys():      
        allPreds['lm'][modTag] = lmPredict(isoPosDict, probsDict[modTag])
    
    #Add collect constitutive mod hits
    for modTag in constDict.keys():
        allPreds['cs'][modTag] = predictProteins(protSeqs, constDict[modTag], eCut, tempDir, cpus)
    
    #Predict relevant proteins
    protHits = {}
    for modTag in protDict.keys():
        protHits[modTag] = predictProteins(protSeqs, protDict[modTag], eCut, tempDir, cpus)
    
    #Output file labelling 
    callFile = '{0}/modCalls.txt'.format(outputFile)
    scoreFile = '{0}/modScores.txt'.format(outputFile)
    protFile = '{0}/protHits.txt'.format(outputFile)
    protScoreFile = '{0}/protScores.txt'.format(outputFile)
    summaryFile = '{0}/summaryFile.txt'.format(outputFile)
    
    
    #Get gaps
    allGaps = getPositions(positionLists)
    
    #Modify tRNA sequences
    modSeqs, modScores, protHitPos, protScores = modifySeqs(allPreds, tRNApos, allCMs, modInfo, nameMap, allGaps, protHits)
    
    columns = ['tRNA'] + allGaps
    
    #Output mod predictions
    writeTSV(callFile, modSeqs, columns)
    #Output mod pred stats
    writeTSV(scoreFile, modScores, columns)
    #Output protein hits
    writeTSV(protFile, protHitPos, columns)
    #Output protein hits
    writeTSV(protScoreFile, protScores, columns)
    
    writeSummaryFile(summaryFile, modSeqs, modScores, protHitPos, protScores)
    
    return callFile, scoreFile, protFile, protScoreFile

def writeSummaryFile(outFile, modSeqs, modScores, protHitPos, protScores):
    """Write a summary output file that is easier to read with a machine"""
    
    with open(outFile, 'w') as outF:
        
        outF.write('tRNA name\tsprinzl position\tprediction\tprobability\tpredicted protein\tprotein bit score\n')
        
        for tRNA in modSeqs.keys():
            for pos in modSeqs[tRNA].keys():
                for mod, score, prot, bit in zip(modSeqs[tRNA][pos].split(','), 
                                                 modScores[tRNA][pos].split(','), 
                                                 protHitPos[tRNA][pos].split(','),
                                                 protScores[tRNA][pos].split(',')):
                    
                    outF.write('tRNA-{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(tRNA, pos, mod, score, prot, bit))
                    
        outF.close()
            
    

if __name__ == '__main__':
    main()
