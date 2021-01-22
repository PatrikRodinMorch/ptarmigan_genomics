#!/usr/bin/python

#%%
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import allel; print('scikit-allel', allel.__version__)
# %%
# Read in Fst output from vcftools into a pandas dataframe as tab delimited
fst = pd.read_csv("/home/prm/Desktop/Lago_small_project/Lago_muta_scot_FST/Fst_15KBwindow.windowed.weir.fst", sep = "\t")

# %%
# Number of non-overlapping genomic windows
len(fst)

#%%
fst['N_VARIANTS'].sum()
# %%
# Convert negative Fst values to 0
fst['MEAN_FST'] = fst['MEAN_FST'].mask(fst['MEAN_FST']<0, 0)
fst['WEIGHTED_FST'] = fst['WEIGHTED_FST'].mask(fst['WEIGHTED_FST']<0, 0)
fst.head(5)

#%%
fst = fst[(fst['CHROM'] != "ChrZ") & (fst['CHROM'] != "ChrW")]

fst.tail()

#%%
max(fst['MEAN_FST'])

# %%
# Convert Fst to the Z score of FST (ZFST)
fst['ZFST'] = (fst['MEAN_FST'] - fst['MEAN_FST'].mean()) / fst['MEAN_FST'].std(ddof = 0)
fst.head(5)

# %%
fst['ZFST'] = fst['ZFST'].mask(fst['ZFST'] < 0, 0)
fst.head(5)


# %%
# Write all ZFST values for manhattan plot
fst.to_csv('Muta_scot_ZFST_all_15kb_autosomes.txt', sep = "\t")

#%%
min(fst['ZFST'])

# %%
# Create a new dataframe with only rows for ZFST>3
ZFST_outliers = fst[fst['ZFST'] >= 3]

#%%
len(ZFST_outliers)

#%%
ZFST_outliers['seqid'] = ZFST_outliers["CHROM"].str.replace("Chr", "")
# %%
# Control the format
ZFST_outliers.head(5)

# %%
# Control that no entries contain any ZFST<3
any(ZFST_outliers['ZFST'] < 3)
 
# %%
len(ZFST_outliers)

# %%
len(ZFST_outliers[ZFST_outliers['ZFST'] >= 6])

# %%
min(ZFST_outliers['ZFST'])

# %%
max(ZFST_outliers['ZFST'])

# %%
## Now we retrieve the annotations for the outlier regions from the Chicken assembly version 6 to which the ptarmigan reads were mapped
annotations = allel.gff3_to_dataframe("Gallus_gallus.GRCg6a.99.chr.gff3.gz", attributes = ['type', 'start', 'end', 'ID', 'biotype', 'description'])

# %%
annotations.head(6)


# %%
len(annotations)

#%%
start_fst = ZFST_outliers['BIN_START'].tolist()
end_fst = ZFST_outliers['BIN_END'].tolist()
seqid_fst = ZFST_outliers['seqid'].tolist()
ZFST = ZFST_outliers['ZFST'].tolist()
chrom_fst = ZFST_outliers['CHROM'].tolist()

# %%
start_anno = annotations.iloc[:, 3].tolist()
end_anno = annotations.iloc[:, 4].tolist()
ID_anno = annotations.iloc[:, 11].tolist()
gene_anno = annotations.iloc[:, -1].tolist()
seqid_anno = annotations.iloc[:, 0].tolist()
type_anno = annotations.iloc[:, 2].tolist()


#%%
### Extracts genes within, and 5Kb above and bellow each window.
test = []

for i in range(len(seqid_anno)):
    for n in range(len(seqid_fst)):
        if (seqid_fst[n] == seqid_anno[i]) & (start_anno[i] >= (start_fst[n] - 15000)) & (end_anno[i] <= (end_fst[n] + 15000)):
            l_fst = [chrom_fst[n], str(start_fst[n]), str(end_fst[n]), str(ZFST[n])]
            l_anno = [seqid_anno[i], str(start_anno[i]), str(end_anno[i]), type_anno[i], ID_anno[i], gene_anno[i]]
            l_comb = l_fst + l_anno
            test.append(l_comb)

#%%
annotation_and_FST = pd.DataFrame(test, columns= ['chrom', 'start_bin', 'end_bin', 'ZFST', 'seqID', 'start_anno', 'end_anno', 'type_anno', 'ID', 'gene'])

#%%
len(annotation_and_FST)


#%%
genes = annotation_and_FST[annotation_and_FST.iloc[:,-1] != '.']
ID_for_GO = annotation_and_FST[annotation_and_FST.iloc[:, 7] == 'gene']


#%%
len(genes)

#%%
len(ID_for_GO)

#%%
genes

#%%
ID_for_GO



# %%
annotation_and_FST.to_csv('Muta_scot_outliers_annotations_15kb_15kboverlap_autosomes.txt', sep = "\t")
genes.to_csv('Muta_scot_outliers_with_genes_15kb_15kboverlap_autosomes.txt', sep = "\t")
ID_for_GO.to_csv('Muta_scot_gene_ID_for_GO_15kb_15kboverlap_autosomes.txt', sep = "\t")


### Create a manhattan plot in matplotlib
# %%
