
import numpy as np
import glob
import pickle
import re
import pandas as pd

'''
    CreateEdgeNetwork - Create Edge file
'''

outputDir = 'data'
f = glob.glob(outputDir + '/*.npz')
assert len(f) > 0, 'gimme files - did not find any %g' % len(f)

for fi in f:
    s = np.load(fi)
    print('==================================\n Processing %s' % fi)
    adjMatrixBootstrap = s['adjMatrixBootstrapQ']
    cost = s['cost']
    assert np.all(adjMatrixBootstrap.shape == cost.shape)
    # we remove significant pairs that are not symmetric
    print('Is adjacency matrix symmetric?', np.allclose(adjMatrixBootstrap, adjMatrixBootstrap.T))
    adjMatrixBootstrap = np.logical_and(adjMatrixBootstrap, adjMatrixBootstrap.T)
    assert np.allclose(adjMatrixBootstrap, adjMatrixBootstrap.T)
    print('Adjacency matrix now symmetric after applying AND transformation')
    G = cost.shape[0]
    nt = G*(G-1)  # number of tests without diagonal
    print('Sparseness %f' % (adjMatrixBootstrap.sum() / float(nt)))
    # Get gene names
    strDir = '/home/mqbssaby/expresults/OscopeBootstrap/data'
    dataGene = None
    if(fi[len(outputDir)+1:len(outputDir)+4] != 'SVZ'):
        datal = ['s123res', 's45res']
        dataList = pickle.load(open(strDir + '/WaterfallData.p', "rb"), encoding='latin1')
        for d in datal:
            if(re.search(d, fi)):
                print('Waterfall data type %s' % d)
                dataGene = dataList[d]
    else:
        print('SVZ data')
        datal = ['qNSC', 'aNSC', 'NB']
        for d in datal:
            if(re.search(d, fi)):
                strFile = strDir + '/SVZprocessedData/SVZ_' + d + '_Rescaled.csv'
                print('Loading file %s' % strFile)
                dataGene = pd.read_csv(strFile, index_col=0)
    assert(dataGene is not None)
    assert(dataGene.shape[0] == G)
    geneNames = dataGene.index.values
    fCreateEdgeRepr = True
    strGeneNameFile = ','.join(map(str, geneNames))  # add commas
    with open(fi[:-4] + '_geneNames.txt', 'w') as f:
        f.write(strGeneNameFile)

    if(fCreateEdgeRepr):
        # Create edge representation
        # G_i, G_j, cost for all significant genes
        nSignificantPairs = adjMatrixBootstrap.sum() / 2.  # symmetric matrix
        assert(nSignificantPairs.is_integer())
        edgeNetwork = []  # np.empty((int(nSignificantPairs), 3), dtype='string, string, float64')
        iterC = 0
        for i in range(G):
            for j in range(i+1, G):
                if(adjMatrixBootstrap[i, j] == 1):
                    #                 edgeNetwork[iterC, :] = [geneNames[i], geneNames[j], cost[i, j]]
                    edgeNetwork.append([geneNames[i], geneNames[j], cost[i, j]])
                    iterC += 1
        a = pd.DataFrame(data=edgeNetwork, columns=['gene1', 'gene2', 'cost'])
        # then rank by cost - distance, the bigger the worse - try to minimize, so order is ascending
        b = a.sort_values('cost')
        strFileNetwork = fi[:-4] + '_edgeNetwork.csv'
        print('Writing network file %s with %g significant gene pairs ' % (strFileNetwork, iterC))
        b.to_csv(strFileNetwork)
        # dump entire cost with gene names
        a = pd.DataFrame(data=cost, columns=geneNames, index=geneNames)
        a.to_csv(fi[:-4] + '_fullCostWithNames.csv')

    # could also add p-values
