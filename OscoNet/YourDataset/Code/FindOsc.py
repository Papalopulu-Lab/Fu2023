import sys
from OscopeBootstrap import qvalue
from OscopeBootstrap.create_edge_network_represention import create_edge_network_representation
from OscopeBootstrap.SyntheticDataset import GetSimISyntheticData, true_adj_matrix
from OscopeBootstrap.oscope_tf import bootstrap_hypothesis_test, get_accuracy, get_metrics_for_different_qvalue_thresholds
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import time

def run_osconet(fileinput,fileoutput,n_bootstrap,verbose=False):


    if (verbose):
        print("Inputfile={0}".format(fileinput))
        print("Outputfile={0}".format(fileoutput))
        #print("n_bootstrap={0}".format(n_bootstrap),flush=True)
        print("n_bootstrap={0}".format(n_bootstrap))

    data_df = pd.read_csv(fileinput)
    genenames = data_df.iloc[:, 0] 
    D = data_df. iloc[:, 1:data_df.shape[1]] 
    D_new =D.rename(index=genenames)
                                               
    grid_points_in_search = 10  # grid size for phase shift parameter estimation., 
    alpha = 0.001  # significance level
    n = D_new.shape[0]
    print("n is {0}".format(n))
    adjacency_matrix, qvalues, cost_matrix, psi_ng = bootstrap_hypothesis_test(n_bootstrap, D_new.values[0:n,:], alpha=alpha,grid_points_in_search=grid_points_in_search)
    print("shape of psi_ng is {0}".format(psi_ng.shape))
    np.savetxt(fileinput[0:-4] + ".psi_ng.csv", psi_ng, delimiter=",")
    print("saving psi_ng to {0}".format(fileinput[0:-4] + ".psi_ng.csv"))
    print(adjacency_matrix.shape)
    print(qvalues.shape)

    G=adjacency_matrix.shape[1]

    gene_names=genenames.values.tolist()
    #add psi_ng extraction into a file!
    edge_network = create_edge_network_representation(adjacency_matrix, 1/cost_matrix, gene_names[0:n])

    edge_network.to_csv(fileoutput, index = False, header=True)

    return edge_network

start_time = time.time()
fileinput = sys.argv[1]
fileoutput=fileinput[0:-4] + '.out.csv'
#fileoutput = sys.argv[1] + '.out.csv'
boot_strap = int(sys.argv[2])
print("CIAOOOOO sboot_strap is {0}".format(boot_strap))
run_osconet(fileinput,fileoutput,boot_strap,verbose=True) 
start_time = time.time()
