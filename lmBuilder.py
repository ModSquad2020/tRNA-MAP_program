#!/usr/bin/env python
# coding: utf-8

# In[23]:


def parseInput():
    """Parse command line input and output files"""
    import argparse
    
    #Argparse setup
    argParser = argparse.ArgumentParser(description = 'This program generates track hubs that display Modomics data')
    argParser.add_argument('-d', '--data_files', required = True, help = 'data reference file')
    argParser.add_argument('-m', '--modification', required = True, help = 'modification short name')
    argParser.add_argument('-p', '--position', required = True, help = 'sprinzl alignment position')
    argParser.add_argument('-o', '--output_file', required = True, help = 'sprinzl alignment position')
    argParser.add_argument('--ignore_positions', required = False, nargs = '*', default = ['74', '75', '76'], help = 'ignore specified positions. Default ignores 74, 75 and 76 (CCA tail) due to random sequences missing it')
    argParser.add_argument('--rsquared_cutoff', required = False, default = 0.15, help = 'R-squared cutoff for individual columns')
    argParser.add_argument('--exclude_species', required = False, nargs = '*', help = 'exclude species, specified by species name in data file. ex: Streptomyces_griseus')
    
    
    kingdoms = {'E': 'Eukaryota', 'A': 'Archaea', 'B': 'Bacteria'}
    
    #Example input: '-m s4U -p 8 --rsquared_cutoff 0.05 --exclude_species Streptomyces_griseus -o ./lm_build_testOutput.txt -d ../data/top5_bact_gold-standard_Modomics.txt'.split()
    clArgs = argParser.parse_args()
    mod = clArgs.modification
    pos = clArgs.position
    data = clArgs.data_files
    igPos = set(clArgs.ignore_positions)
    lmCutoff = float(clArgs.rsquared_cutoff)
    outFile = clArgs.output_file
    exclSp = clArgs.exclude_species
    
    return mod, pos, data, lmCutoff, outFile, igPos, exclSp

def parseData(inputFile):
    """Parse Modomics data tsv format file"""
    
    data = {} #Store data as a dictionary of {seqName: [sequence]}
    
    listKey = [] #Store positioning of each part as a list
    
    with open(inputFile, 'r') as inF:
        
        lines = inF.readlines()
        
        #Header and tsv column info
        header = lines[0].strip().split('\t')
        listKey = header[4:]
        
        for line in lines[1:]:
            
            splitLine = line.strip().split('\t')
            
            #Name sequence
            seqName = '-'.join(splitLine[0:2])
            
            seq = splitLine[4:]
            
            data[seqName] = seq
            
    return data, listKey

def cleanNfilter(data, posKey, mod, refPos, exclSp):
    """Removes extra modifications and filters out unmodifiable nucleotides"""
    #Import from parent dir since I'm running in a subdirectory
    import sys
    sys.path.append('..')
    
    from tRNAinfo import short2original as sto
    
    #Store cleaned and filtered data
    cfDict = {} #Array of data
    colArrays = [] #Array of data in each column of alignment
    truthArray = [] #Array of truth
    nameArray = [] #Array of names ordered the same as in truth
    
    refBase = sto[mod]
    
    #Iterate through each sequence in the dataset
    for seqName in data.keys():
        
        newSeq = []
        truth = 0
        
        #Skip over specified species
        sp = seqName.split('-')[0]
        
        if sp in exclSp:
            pass
        else:
        
            #Iterate through each position and back-translate to original base
            for pos, resi in zip(posKey, data[seqName]):

                # Update the truth if the correct mod is present.
                if resi == mod and pos == refPos:
                    truth = 1

                #Handle gaps. Really would work better if I added them to the dictionary. But feeling lazy lol   
                try:
                    newSeq.append(sto[resi])
                except KeyError:
                    if resi == '-':
                        newSeq.append(resi)
                    else:
                        pass

                #If we are on the reference position, and the correct unmodified base is present, 
                # add seq name to the cleaned and filtered dictionary

                if resi == '-': #Again, not including the - is really screwing me!
                    pass
                elif pos == refPos and sto[resi] == refBase:
                    cfDict[seqName] = []
                else:
                    pass
        
            try:
                cfDict[seqName]+=newSeq
                truthArray.append(truth)
                nameArray.append(seqName)
            except KeyError:
                pass
    
    
    #Build column-wise datasets
    colArrays = {x:[] for x in posKey} #Array of data in each column of alignment
    for seq in cfDict.keys():
        
        for pos, resi in zip(posKey, cfDict[seq]):
            
            colArrays[pos].append(resi)
            
    return cfDict, truthArray, nameArray, colArrays
    
def writeLM(outFile, modelColNames, modelCols, Truth, names, mod, pos):
    """Write tsv output file for linear model data"""
    
    with open(outFile, 'w') as outF:
        
        #Write header
        outF.write( 'seqName\t' + mod +'_'+ pos + '\t' + '\t'.join(modelColNames) + '\n')
        
        for name, tru, colVals in zip(names, Truth, modelCols):
            outF.write ( name + '\t' + str(tru) + '\t' + '\t'.join(colVals) + '\n')
            
        outF.close()
    
def lmBuild(lmData, Truth, names, cols, R2cut, colKey, skipCols): 
    """Builds linear model from linear-regression training data; searches tRNAs against linear model"""
    import numpy as np
    from sklearn.preprocessing import OneHotEncoder as OHE
    from sklearn.linear_model import LinearRegression as LR
    
    #Add prediction set to predictor set for encoding.
    #Iterate through and predict tRNAs
    
    predVars = np.array(lmData)
    predState = np.array(Truth)
    
    colPreds = []
    
    #Remove specified columns
    colKeys = set(cols.keys())
    colKeys = colKeys - skipCols
        
    scoreDict = {} #Store prediction scores for each column (score_value: [columns that have score])
                   #Just so its easy to sort!
    
    #First: Build individual linear models for each predictor variable
    for pos in colKeys:
        
        #Convert to vertical vector
        colCol = [ [x] for x in cols[pos]]
        predCol = np.array(colCol)
        
        #Encode prediction column in onehot
        encoding = OHE(sparse=False)
        encodedArray = encoding.fit_transform(predCol)
        
        #Instantiate linear model  ###########################NOTE: I MAY HAVE THIS BACKWARDS. NOT EXACTLY SURE WHICH WAY DATA GETS FED INTO LM!!!
        linearModel = LR().fit(encodedArray, predState)
        
        #Return the R^2 of the prediction
        try:
            scoreDict[linearModel.score(encodedArray, predState)].append(pos)
        except KeyError:
            scoreDict[linearModel.score(encodedArray, predState)] = [pos]
            
        #print(linearModel.coef_)
        #print(linearModel.intercept_)
        
    #Add together all columns that make the cutoff
    lmPos = []
    for score in sorted(scoreDict.keys(), reverse = True):
        
        #Add position of above cutoff
        if score >= R2cut:
            lmPos += scoreDict[score]      
        else:
            pass
    
    newModel = []
    orderedLMpos = []
    #Iterate through each sequence, make sure lm positions are in correct spot
    for seqName in lmData.keys():
        
        #Store relevant sequence positions
        newSeq = []
        singleLMpos = [] #Store lm columns in correct order!
        for pos, base in zip(colKey, lmData[seqName]):
            
            #Add nucleotide to new sequence if in predictive column
            if pos in lmPos:
                newSeq.append(base)
                singleLMpos.append(pos)
            else:
                pass
                
        newModel.append(newSeq)
        orderedLMpos.append(singleLMpos)
    
    #New data array
    newArray = np.array(newModel)
    
    #Encode in one hot
    encoding = OHE(sparse=False)
    newEncode = encoding.fit_transform(newArray)
    
    #Build new linear model
    combinedLM = LR().fit(newEncode, predState)
    
    #Score new linear model:
    combinedScore = combinedLM.score(newEncode, predState)
    
    return combinedScore, newModel, orderedLMpos[0]
        
    
def main():
    """ Execute functions to identify best linear model """
    
    mod, pos, dataFile, cut, outputFile, ignores, excluded = parseInput() #add in ignores
    
    # Parse input file
    data, posKey = parseData(dataFile)
    
    # Filter out appropriate sequences
    niceData, truth, names, cols = cleanNfilter(data, posKey, mod, pos, excluded)
    
    # Build linear regression
    score, modelCols, modelColNames = lmBuild(niceData, truth, names, cols, cut, posKey, ignores) 
    
    writeLM(outputFile, modelColNames, modelCols, truth, names, mod, pos)

if __name__ == '__main__':
    main()


# In[ ]:




