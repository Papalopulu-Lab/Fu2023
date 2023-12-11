# Network 2022

The scripts in this sub-directory can be used to build a network of transcriptional interactions with miRNA regulation, using the latest databases available as of May 2023. 

The final constructed network can be readily be provided by the authors upon request. The identified regulatory motifs are provided in this repository. The dynamic_network_motifs.txt file contains the genes predicted to oscillate in this study.

# Requirements

We recommend to use an anaconda distribtuion of python3 as this will come with all the requisite libraries pre-installed.

# Instructions

## Download Data
### ReMap Data

Download ReMap2022 h38 All peaks and CRM data (https://remap2022.univ-amu.fr/) into your Network Construction directory.

### Ensembl Genes Data

Download Human Ensembl Genes data from Ensembl BioMart (https://www.ensembl.org/biomart/martview/) as biomart.txt into your Network Construction directory. Required fields: 

- Gene stable ID
- Gene stable ID version
- Transcript stable ID
- Transcript stable ID version
- Exon stable ID
- Chromosome/scaffold name
- Gene start (bp)
- Gene end (bp)
- Strand, Transcript start (bp)
- Transcript end (bp)
- Transcription start site (TSS)
- Transcript length (including UTRs and CDS)
- Gene name
- Transcript name
- Gene type
- Transcript type
- 5' UTR start
- 5' UTR end
- 3' UTR start
- 3' UTR end
- Exon region start (bp)
- Exon region end (bp)
- Exon rank in transcript

### Ensembl miRNA Data

Download Human miRNA data from Ensembl BioMart (https://www.ensembl.org/biomart/martview/) as ensembl_miRs.txt into your Network Construction directory. Filter for Gene Type = miRNA. Required fields: 

- Gene stable ID
- Gene stable ID version
- Transcript stable ID
- Transcript stable ID version
- Unspliced (transcript) sequence 

### miRBase Data

Download hsa.gff3, hairpin.fa and mature.fa files Release 22.1 from miRBase (https://www.mirbase.org/download/) into your Network Construction directory.

### miRTarBase Data

Download hsa_MTI.xlsx Release 8.0 from mirTarBase (https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/php/download.php) into your Network Construction directory.

## Control File

Some scripts work by reading from a control file. A copy of the control file can also be found in this repository. Ensure you have edited the paths to match those on your system and saved the file in your Network Construction directory.

```
#file_locations_for_map_construction

@	INFILE_DIRECTORY=/home/NetworkConstruction_2022/

@	MIRTAR_FILE=/home/NetworkConstruction_2022/hsa_MTI.txt

@	MIRHOSTS=hsa_hosts.txt

@	GENOME_FILE=biomart.txt

@	REMAP_FILE=/home/NetworkConstruction_2022/remap2022_all_macs2_hg38_v1_0-REFORMATED-crm_filtered.bed

@	NETWORK_DIRECTORY=/home/NetworkConstruction_2022/

@	MAP_FOLDER=/home/NetworkConstruction_2022/work_files/map/

@	MIRTAR_edge_file=/home/NetworkConstruction_2022/mirtar_edges.txt
```

## Data 
### Reformat ReMap Data 

Before using the ReMap TFBS data, reformat it using REMAP2022.py.

```
python REMAP2022.py [remap2022_all_macs2_hg38_v1_0.bed]

Output: remap2022_all_macs2_hg38_v1_0-REFORMATTED.bed
```

### CRM Filtering

ReMap TFBS data should be ﬁltered to binding sites which are located in the CRMs, to remove non-speciﬁc binding. Perform this step by running the python script in chip_in_crm in the console/terminal with the CRM and All Peaks files as positional arguments. 

```
python chip_in_crms.py [remap2022_crm_macs2_hg38_v1_0.bed] [remap2022_all_macs2_hg38_v1_0-REFORMATTED.bed]

Output: remap2022_all_macs2_hg38_v1_0-REFORMATTED-crm_filtered.bed
```

## Assigning Transcriptional Regulator Binding Site (TRBS) to target genes

ReMap TFBS data should be assigned to the Ensembl TSS which is closest within a 60kb region (10kb downstream and 50kb upstream of TSSs) using peak_MULTI.py. Where two or more TSSs are within this range, the TSS closest will be assigned to that TFBS. The upstream and downstream distance can be edited on line 100 of the script. Distances between each TSS and TFBS are recorded. Run this script using the Control File as a positional argument. Setting ulimit -n to 200000 may improve performance.

    python peakMULTI.py [control_file.txt] ulimit -n 200000

    Output: /mirHOST/ & /work_files/genome_split/

## Assigning microRNAs to host transcripts
### Identifying microRNA host genes

As microRNAs can be produced at the same time as host genes (Morlando et al., 2008) microRNAs are assigned regulatory inputs from the transcripts with which their annotations overlap. As the microRNA coordinates will be for the mature microRNA, this will also allow the mature transcripts to gain the regulatory input from their microRNA genes. 

Overlap with transcripts and positions within the transcripts can be performed using microME_2022.py. This script will
access the microRNA loci as annotated in miRBase and look for overlapping annotations between the microRNA and genes in the biomart download.

This script takes three positional arguments, (1) the species code "hsa" for Homo sapiens, (2) the Biomart ﬁle containing the transcript positions and (3) the miRBase hsa.gff3 file. 

```
python microME_plus2022.py [hsa] [biomart.txt] [hsa.gff3]

Output: /mirHOST/hsa_boxdata.txt & /mirHOST/hsa_hosts.txt & mirHOST/hsa_structure.txt
```

### Assigning host gene regulation to microRNAs

MicroRNAs next need to be assigned regulatory input from their host genes and the connections incorporated into the network. MicroRNAs within the peak_MULTI output ﬁles will then inherit a copy of the TSS regulation of the genes in which they are contained, while also maintaining the TSS as listed in the BioMart data.

mir_chip_MULTI.py takes ﬁve positional arguments. (1) The ccm directory, (2) the microRNA host ﬁle produced in the previous step, (3) the biomart ﬁle, (4) a ﬁle name to write TF ids to (this is a writeable path) and (5) the reformatted ReMap ﬁle. This script also requires a ﬁle called "unspec_chip.txt" to be created in the work_ﬁle directory. This ﬁle should contain the names of any factors which should not be included in the analysis. An example ﬁle is included with the scripts. Factors were excluded based on the non-speciﬁcity of the antibodies used, which would identify multiple genes. For example "RUNX" which identiﬁes binding sites for both RUNX1 and RUNX2.

    python mir_chip_multi.py [work_files/chip_close_match] [mirHOST/hsa_hosts.txt] [biomart.txt] [home/NetworkConstruction_2022/tf_id_file.txt] [remap2022_all_macs2_hg38_v1_0-REFORMATED-crm_filtered.bed]

    Output: tf_id_file.txt
    
## Generating the transcriptional and translational network
### Reformatting miRTarBase data

miRTarBase data should be converted from an excel ﬁle to a tab-delimited ﬁle using excel (save as). This ﬁle will then be reformatted into a structure matching the networks using miRTarToMap.py. This script requires two positional arguments, (1) the biomart_ﬁle path and (2) the microRNA host ﬁle path produce earlier.

```
python miRTarToMap.py [biomtart.txt] [hsa_MTI.txt]

Output: mirtar_edges.txt
```

### Collecting all the map files together

This step collects all of the separate map ﬁles in the map folder and joins them together with the miRTarBase ﬁle to produce a connected network.

NOTE: The network is not complete at this step and requires the correction of the miRTarBase codes in miRBase codes.

Joining of the network can be performed using build_map_2022.py. This script takes the Control File as a positional argument.

    python build_map2022.py [control_file.txt]

    Output: REMAP2022_map.txt
    
### Analysing peaks by experiment 

The next two scripts collapse the network, removing unecessary information. The network is simpliﬁed to reduce to the multiple connections between two nodes to an edge weight. This first script takes, one positional argument which is the path to the output from the previous step "REMAP2022_map.txt". All duplicate edges are collapsed and counted to make an edge weight. The distance to the TSS reported in this version of the network is the smallest distance of any TR binding site for that edge.

```
python analyse_peaks_by_experiment.py [REMAP2022_map.txt]

Output: REMAP2022_select-collapsed.txt
```

### Collapsing select to mean

This script then takes the output from above as a positional argument to complete the collapsing process.

```
python collapse_select_to_mean.py [REMAP2022_select-collapsed.txt]
    
Output: REMAP2022_select-collapsed_ctm.txt
```

## Filtering Network

For each dataset, the number of binding sites for each TR will be recorded at each TSS. Following the removal of outliers (> 2 SD of the mean), the mean number of binding sites for each TR is calculated. This is performed at the level of the individual ChIP experiments to avoid the inﬂuence of variation between datasets. If the average number of binding sites at TSSs for a given TR is n, then any TSS where fewer than n sites are observed will be removed from the network. Means for interactions are calculated after removing extreme outliers. TSS, where ≥n sites are found, will be maintained. This should ﬁlter out any non-speciﬁc TFBS interactions leaving only high conﬁdence sites remaining.

The first step of filtering is performed using edge_weight_analysis.py on the collapsed network generated in the previous step. 

```
python edge_weight_analysis.py [REMAP2022_select-collapsed_ctm.txt]

Output: REMAP2022_map_select-collapsed_ctm-edge_weight_stats80th.txt
```

This generates an Edge Weight Stats file which can be used to run the 2nd part of filtering.

    python edge_weight_filter.py [REMAP2022_select-collapsed_ctm.txt] [REMAP2022_map_select-collapsed_ctm-edge_weight_stats80th.txt]

    Output: REMAP2022_map_select-collapsed_ctmEdgeFmean2SDoM.txt

    
## Correcting miRNA codes 
### Adding miRBase gene codes

The current network contains miRTarBase interaction codes in place of a real gene ID for the microRNA. To ﬁx this, the network is scanned for microRNA, and the gene codes are replaced for miRBase gene codes. This is performed using correct_mir_codes.py. This script takes two positional arguments (1) the REMAP2022_map.txt ﬁle path, (2) the miRBase mature.fa ﬁle.

```
python correct_mir_codes.py [REMAP2022_map_select-collapsed_ctmEdgeFmean2SDoM.txt] [mature.fa]

Output: REMAP2022_map_select-collapsed_ctmEdgeFmean2SDoM-fixed_mirs.txt
```

### Restoring header 

This is merely a cosmetic script, to return the header titles to the network to allow for easier navigation of the completed network.

    python headerer.py [REMAP2022_map_select-collapsed_ctmEdgeFmean2SDoM-fixed_mirs.txt]

    Output: REMAP2022_map_select-collapsed_ctmEdgeFmean2SDoM-fixed_mirs_header.txt
    
### Converting miRBase to ensembl codes 

Now both the miRBase and miRTarBase data have been given miRBase codes. The final step is to convert these to Ensembl codes so they match up with the ReMap miRNA data. This way we can look for feedback interactions during the motif analysis. miRBase mature miRNA codes are first converted from to hairpin miRNA codes, as Ensembl and thus ReMap miRNA data consists of hairpin miRNA entries. Then hairpin sequences from miRBase are cross-referenced with the ensembl_miRs file to build a dictionary assigning Ensembl codes to each hairpin miRNA in miRBase. Then this dictionary is used to parse through the network and convert all miRNAs to Ensembl codes.

This script takes 4 positional arguments (1) the mature.fa and (2) hairpin.fa files from miRBase, (3) the ensembl_miRs.txt file and (4) the network file.

```
python miRBase-ensembl_converter.py [mature.fa] [hairpin.fa] [ensembl_miRs.txt] [REMAP2022_map_select-collapsed_ctmEdgeFmean2SDoM-fixed_mirs_header.txt]

Output: REMAP2022_map_select-collapsed_ctmEdgeFmean2SDoM-fixed_mirs_header_converted.txt
```

## Motif Identification
### Reducing Network Size

As all motifs that we are interested in contain feedback, the network can be ﬁltered to contain only target genes which are also regulators. This increases the speed of network motif identiﬁcation. Networks were ﬁltered to contain only interacting components using "just_connections_map.py".

This takes one positional argument, which is the ﬁltered network from the previous step and
produces a new network ﬁle with "just_interactions" added to the ﬁle name.

    python just_connections_map.py [REMAP2022_map_select-collapsed_ctmEdgeFmean2SDoM-fixed_mirs_header_converted.txt]

    Output: REMAP2022_map_select-collapsed_ctmEdgeFmean2SDoM-fixed_mirs_header_converted-just_interactions.txt
    
### Network Motifs

Network motifs are identiﬁed by looking for patterns within the transcriptional-translational network edge list utilising get_motifs_quicker.py. Autoregulatory loops are calculated ﬁrst as other motifs are dependent on these loops. TRs in autoregulatory loops are then used to shortlist the search for further motifs. Motif discovery is based on boolean logic looking for distinct patterns on interaction. For example, for autoregulation instances are identiﬁed where gene-A targets gene-A. Autoregulation with microRNA feedback would be gene-A targets gene-A, gene-A targets microRNA-x and microRNA-x targets gene-A.

The python script get_motifs_quicker.py takes one positional argument which is the network ﬁle path. This script will create a new directory in the same folder as the network called /fast_motifs. In the fast_motifs directory will be individual text ﬁles which contain all of the motifs and a motif summary ﬁle containing all of the numbers of instances of each motif within the network.

```
python get_motifs_quicker.py [REMAP2022_map_select-collapsed_ctmEdgeFmean2SDoM-fixed_mirs_header_converted-just_interactions.txt]

Output: /fast_motifs/
```

## Network 2022 Complete

Congratulations! You have built the Network. You're now ready to mine for interesting dynamic molecules, oscillatory transcriptional factors, regulatory feedbacks motifs and much more.
