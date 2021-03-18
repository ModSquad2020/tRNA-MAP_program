"""
Use Bokeh to visualize mod predictions

input files:
- tRNA mod calls tab-delimited text file
- tRNA mod scores text file
- protein prediction file

outputs:
- interactive html file bilt using bokeh

To Add:
- matt wants a confidence measure for high bit scores...
- Add in protein info!
- Add in color key
- 3' and 5' ends!
"""

def parseInput():
    """Parse command line input and output files"""
    import argparse
    
    #Argparse setup
    argParser = argparse.ArgumentParser(description = 'This program generates track hubs that display Modomics data')
    argParser.add_argument('-c', '--mod_calls', required = True, help = 'modification calls file')
    argParser.add_argument('-o', '--output_file', required = True, help = 'output html file')
    #argParser.add_argument('-k', '--kingdom', required = False, choices = {'E', 'A', 'B'}, help = 'Domain of organism')
    argParser.add_argument('-s', '--mod_scores', required = True, help = 'modification scores file')
    argParser.add_argument('-e', '--enzyme_predictions', required = True, help = 'enzyme predictions file')
    
    kingdoms = {'E': 'Eukaryota', 'A': 'Archaea', 'B': 'Bacteria'}
    #'-c ./strepneumo_test/modCalls.txt -s ./strepneumo_test/modScores.txt -e ./strepneumo_test/output.tsv -o ./strepneumo_test'.split()
    clArgs = argParser.parse_args()
    callFile = clArgs.mod_calls
    #orgKing = kingdoms[clArgs.kingdom]
    outFile = clArgs.output_file
    scoreFile = clArgs.mod_scores
    protFile = clArgs.enzyme_predictions
    
    return callFile, scoreFile, protFile, outFile

def readTSV(inputFile, header = True):
    """Read a TSV with row and column headers"""
    
    #Store mod info: mod: {modified_position: {'+': [tRNAs], '-': [tRNAs]]}
    seqsDict = {}
    
    with open(inputFile, 'r') as inF:
        
        #Handle header presence
        rawLines = inF.readlines()
        
        #TSV columns
        cols = []
        lines = []
        if header == True:
            cols = rawLines[0].strip().split('\t')
            lines = rawLines[1:]
        else:
            lines = rawLines
        
        #parse each line
        for line in lines[1:]:
            splitLine = line.strip().split('\t')
            
            #tRNA name
            tRNAname = splitLine[0]
            
            #Make dictionary of sprinzl positions
            seqsDict[tRNAname] = {}
            
            #Iterate through sprinzl positions
            for sprinzlPos, base in zip(cols[1:], splitLine[1:]):
                
                seqsDict[tRNAname][sprinzlPos] = base
                    
        inF.close()
    
    return seqsDict

def parseProts(inputFile):
    """Parse protein file, return dictionary sorted by modification and position"""
    
    hitDict = {}
    
    with open(inputFile) as inF:
        
        lines = inF.readlines()
        
        for line in lines[:-2]:
            
            #Handle entries that dont have associated mod and position
            try:
                splitLine = line.strip().split('\t')
                mod = splitLine[3]
                pos = splitLine[4].split(' ')
                domains = splitLine[0].split(',')
                accNo = splitLine[1].split(',')
                
                eVals = splitLine[2].split(',')
                
                for p in pos:
                    
                    for d, e, acc in zip(domains, eVals, accNo):
                        
                        try:
                            try:
                                hitDict[p][mod][d] = [e, acc]
                                    

                            except KeyError:
                                hitDict[p][mod] = {d:[e, acc]}

                        except KeyError:
                            hitDict[p] = {mod: {d:[e, acc]}}
                    
            except IndexError:
                pass
                   
    return hitDict


def bit2prob(bitScore):
    """Convert bit score to probability, assuming a binary outcome"""
    
    return 1/(1+2**(-bitScore))
    
def hexCol(prob):
    """Returns a hexadecimal color that represents """
    return "#%02x%02x%02x" % (235, int(round(235-(235*(prob)/3))), int(round(235-(235*(prob)/3))))

def makeFig(callDict, scoreDict, protDict, outF):
    """Make an interactive cloverleaf output figure"""
    import numpy as np
    from bokeh.plotting import figure, show, save
    from bokeh.models import ColumnDataSource, Grid, LinearAxis, Plot, Text, ColorBar, HoverTool
    from tRNAinfo import cloverCoords
    
    """
    Goals for this code
    - Get rid of axes and grid lines
    - adjust sizes of circles for variable loop
    - add key for heat map
    - general fine-tuning
    """
    
    
    #Construct color key for heat map
    colKey = []
    keyX = []
    keyYtop = []
    keyYbottom = []
    keyLabels = [] #ADD LABELS ONCE COLOR STUFF FIGURED OUT!
    
    xPos = 18.539
    #Construct blues
    #for prob in np.arange(0, 0.5, 0.01):
    #    keyYtop.append(-110)
    #    keyYbottom.append(-130)
    #    keyX.append(18 + prob*2/0.01)
    #    colKey.append("#%02x%02x%02x" % (int(round(255-(225*(0.5-prob)))), int(round(255-(160*(0.5-prob)))), 255))
            
    #construct reds
    for prob in np.arange(0, 1, 0.01):
        keyYbottom.append(-110)
        keyYtop.append(-130)
        keyX.append(18 + prob*2/0.01)
        colKey.append(hexCol(prob))
    
    #Construct text for heat map key
    mapTextX = [min(keyX)+((max(keyX)-min(keyX))/2)-1, min(keyX)-1, min(keyX)+((max(keyX)-min(keyX))/2)-1, max(keyX)-1]
    mapTextY = [-100, -145, -145, -145]
    labelText = ['Probability', '0', '0.5', '1']
    
    #Iterate through tRNA 120.418
    for tRNA in callDict.keys():
        
        ################################################################################################################
        #Format data vectors
        
        #Store data vetors
        X = []
        Y = []
        col = []
        scores = []
        bases = []
        probs = []
        positions = []
        sizes = []
        proteins = []
        
        #Add different values to different vectors
        for pos in cloverCoords.keys():
            
            positions.append(pos)
            bases.append(callDict[tRNA][pos])
            
            #Add protein information
            try:
                #generate labels
                protLabels = ['{0}: {1}, {2}'.format(x, protDict[pos][callDict[tRNA][pos]][x][1], 
                                                     protDict[pos][callDict[tRNA][pos]][x][0]) for x in protDict[pos][callDict[tRNA][pos]] ]
                
                proteins.append('; '.join(list(protLabels)))
            except KeyError:
                proteins.append('None')
                                
            
            X.append(cloverCoords[pos][0])
            Y.append(-cloverCoords[pos][1])
            
            modPos = '_'.join([callDict[tRNA][pos], pos])
            
            #Select color
            if scoreDict[tRNA][pos] != '-':
                
                #append probabilities
                prob = bit2prob(float(scoreDict[tRNA][pos]))
                probs.append(prob)
                
                #Handle -inf values which I haven't fixed with laplace pseudocount
                try:
                    scores.append(round(float(scoreDict[tRNA][pos]), 1))
                except OverflowError:
                    scores.append(float(scoreDict[tRNA][pos]))
                    
                #Color the heat map based on mod score
                col.append(hexCol(prob))
                
                #if prob > 0.5:
                #    col.append("#%02x%02x%02x" % (255, int(round(255-(225*(prob-0.5)))), int(round(255-(160*(prob-0.5))))))
                
                #elif prob < 0.5:
                #    col.append("#%02x%02x%02x" % (int(round(255-(225*(0.5-prob)))), int(round(255-(160*(0.5-prob)))), 255))
                
                #else:
                #    col.append("#%02x%02x%02x" % (255, 255, 255))
                    
            else:
                
                probs.append('None')
                scores.append('None')
                col.append("#%02x%02x%02x" % (235, 235, 235))
                
            #Add size of circle
            if pos[0] == 'e':
                sizes.append(cloverCoords[pos][2]+12)
            else:
                sizes.append(cloverCoords[pos][2]+13)
        
        ##########################################################################################################
        
        cloverSource = ColumnDataSource({'x': X, 'y': Y, 'col': col, 
                                         'bases': bases, 'scores': scores, 
                                         'probs': probs, 'positions': positions, 
                                         'sizes': sizes, 'proteins': proteins})
        
        Tooltips = [('modification', '@bases'), 
                    ('score', '@scores'),
                    ('position', '@positions'),
                    ('possible enzymes', '@proteins')]
        
        #Instantiate cloverleaf structure 
        clover = figure(title = '{0}'.format(tRNA), x_axis_label = '', y_axis_label = '',
                        toolbar_location=None, tools = '')#tooltips = Tooltips
        
        #This works with the mouse-over pop ups
        sprinzlCircles = clover.circle('x', 'y', size = 'sizes', fill_color = 'col',
                                          line_color = 'grey', name = 'positions', source = cloverSource)
        
        circleHover = HoverTool(renderers = [sprinzlCircles], tooltips = Tooltips)
        clover.add_tools(circleHover)
        
        #Plot text
        baseText = Text(x = 'x', y = 'y', 
                        text = 'bases', text_color = 'black', 
                        text_align = 'center', text_baseline = 'middle')
        
        clover.add_glyph(cloverSource, baseText)
        
        #color key map
        keyMap = ColumnDataSource({'x': keyX, 'yTop': keyYtop, 'yBottom': keyYbottom, 'col': colKey})
        
        #create color key:
        colorMap = clover.vbar('x', width = 2, top = 'yTop', bottom = 'yBottom', line_alpha = 0, fill_color = 'col', source = keyMap)
        
        labelMap = ColumnDataSource({'x': mapTextX, 'y': mapTextY, 'text': labelText})
        mapText = Text(x = 'x', y = 'y', text = 'text', text_color = 'black', 
                        text_align = 'center', text_baseline = 'middle')
        clover.add_glyph(labelMap, mapText)
        
        #Remove excess visual info from default plotting
        clover.xaxis.major_tick_line_color = None
        clover.yaxis.major_tick_line_color = None
        clover.xaxis.minor_tick_line_color = None
        clover.yaxis.minor_tick_line_color = None
        clover.xgrid.grid_line_color = None
        clover.ygrid.grid_line_color = None
        clover.xaxis.major_label_text_font_size = '0pt'
        clover.yaxis.major_label_text_font_size = '0pt'
            
        save(clover, filename = '{0}/{1}.html'.format(outF, tRNA))
        
def main(inCL = True, calls = None, scores = None, proteins = None, outFile = None):
    import pandas as pd
    from tRNAinfo import sprinzl2coords
    
    #use command line
    if inCL:
        #Recieve input commands
        calls, scores, proteins, outFile = parseInput()
    
    callsDict = readTSV(calls)
    scoresDict = readTSV(scores)
    domainHits = parseProts(proteins)
    
    makeFig(callsDict, scoresDict, domainHits, outFile)