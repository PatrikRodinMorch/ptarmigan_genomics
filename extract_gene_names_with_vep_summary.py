#!/usr/bin/python

#%%
import pandas as pd
import numpy as np
import glob
import os


#%%
def vep_summary(df):
    path = "/home/prm/Desktop/Ptarmigan_analysis/vep_stuff/final/" + df
    df = pd.read_table(path, sep="\t", skiprows=29)
    individual = df['IND'][0]
    delet = len(df[df['SIFT'].str.startswith('deleterious(')])
    lof = len(df[(df['IMPACT'] == 'HIGH') & (df['Consequence'] == 'stop_gained')])
    missense = len(df[df['Consequence'] == 'missense_variant'])
    syn = len(df[df['Consequence'] == 'synonymous_variant'])

    return individual, delet, lof, missense, syn
    
def extract_geneID(df):
    path = "/home/prm/Desktop/Ptarmigan_analysis/vep_stuff/final/" + df
    df = pd.read_table(path, sep="\t", skiprows=29)  
    gene_id_pop_missense = (df[df['Consequence'] == 'missense_variant']['Gene']).unique().tolist()
    gene_id_pop_delet = (df[df['SIFT'].str.startswith('deleterious(')]['Gene']).unique().tolist()
    gene_id_pop_lof = (df[(df['IMPACT'] == 'HIGH') & (df['Consequence'] == 'stop_gained')]['Gene']).unique().tolist()
    
    return gene_id_pop_missense, gene_id_pop_delet, gene_id_pop_lof

def extract_mutation_location(df):
    path = "/home/prm/Desktop/Ptarmigan_analysis/vep_stuff/final/" + df
    df = pd.read_table(path, sep="\t", skiprows=29)   
    individual_miss = df[df['Consequence'] == 'missense_variant']['IND'].tolist()
    individual_del = df[df['SIFT'].str.startswith('deleterious(')]['IND'].tolist()
    individual_lof = df[(df['IMPACT'] == 'HIGH') & (df['Consequence'] == 'stop_gained')]['IND'].to_list()
    location_miss = df[df['Consequence'] == 'missense_variant']['Location']
    location_del = df[df['SIFT'].str.startswith('deleterious(')]['Location']
    location_lof = df[(df['IMPACT'] == 'HIGH') & (df['Consequence'] == 'stop_gained')]['Location']
    
    return individual_miss, individual_del, individual_lof, location_miss, location_del, location_lof

def append_mutations(file, poplabel):
    individual, delet, lof, missense, syn = vep_summary(file)
    pop.append(poplabel)
    ind.append(individual)
    syn_mut.append(syn)
    miss_mut.append(missense)
    delet_mut.append(delet)
    lof_mut.append(lof)
    miss_syn.append(float(missense)/float(syn))
    delet_syn.append(float(delet)/float(syn))
    lof_syn.append(float(lof)/float(syn))

def append_geneID(file, poplabel):
    gene_id_pop_missense, gene_id_pop_delet, gene_id_pop_lof = extract_geneID(file)

    for l in gene_id_pop_missense:
        pop_geneID_miss.append(poplabel)
        miss_id_pop.append(l)

    for l in gene_id_pop_delet:
        pop_geneID_del.append(poplabel)
        del_id_pop.append(l)
     
    for l in gene_id_pop_lof:
        pop_geneID_lof.append(poplabel)
        lof_id_pop.append(l)

def append_location(file, poplabel):
    individual_miss, individual_del, individual_lof, location_miss, location_del, location_lof = extract_mutation_location(file)
    ind_miss = individual_miss
    ind_del = individual_del
    ind_lof = individual_lof

    miss = location_miss.str.split(':', expand = True)
    de = location_del.str.split(':', expand = True)
    lof = location_lof.str.split(':', expand = True)
    miss_location = list(miss.iloc[:,1])
    miss_chrom = list(miss.iloc[:,0])
    del_location = list(de.iloc[:,1])
    del_chrom = list(de.iloc[:,0])
    lof_location = list(lof.iloc[:,1])
    lof_chrom = list(lof.iloc[:,0])

    for i,n,z in zip(ind_miss, miss_chrom, miss_location):      
        id_miss.append(i)  
        chr_miss.append(n)        
        loc_miss.append(z)
        pop_miss.append(poplabel)

    for i,n,z in zip(ind_del, del_chrom, del_location):
        id_del.append(i)  
        chr_del.append(n)        
        loc_del.append(z)
        pop_del.append(poplabel)

    for i,n,z in zip(ind_lof, lof_chrom, lof_location):
        id_lof.append(i)  
        chr_lof.append(n)        
        loc_lof.append(z)
        pop_lof.append(poplabel)

def save_geneID(df, pop):
    gene_id_pop_missense, gene_id_pop_delet, gene_id_pop_lof = extract_geneID(df)
    geneID_pop_columns = ['GeneID']
    geneID_pop_missense = pd.DataFrame(gene_id_pop_missense, columns = geneID_pop_columns)
    #geneID_pop_missense_uniq = pd.DataFrame(geneID_pop_missense['GeneID'].unique())
    geneID_pop_missense.to_csv(pop+'_geneID_pop_missense_uniq_chrZ.txt', sep = '\t', header = True, index = False)
   
    geneID_pop_delet = pd.DataFrame(gene_id_pop_delet, columns = geneID_pop_columns)
    #geneID_pop_delet_uniq = pd.DataFrame(geneID_pop_delet['GeneID'].unique())
    geneID_pop_delet.to_csv(pop+'_geneID_pop_delet_uniq_chrZ.txt', sep = '\t', header = True, index = False)
    
    geneID_pop_lof = pd.DataFrame(gene_id_pop_lof, columns = geneID_pop_columns)
    #geneID_pop_lof_uniq = pd.DataFrame(geneID_pop_lof['GeneID'].unique())
    geneID_pop_lof.to_csv(pop+'_geneID_pop_lof_uniq_chrZ.txt', sep = '\t', header = True, index = False)

#%%
ind = []
pop = []
syn_mut = []
miss_mut = []
delet_mut = []
lof_mut = []
miss_syn = []
delet_syn = []
lof_syn = []
pop_geneID = []
pop_geneID_miss = []
pop_geneID_del = []
pop_geneID_lof = []
miss_id_pop = []
del_id_pop = []
lof_id_pop = []
loc_miss = []
chr_miss = []
id_miss = []
loc_del = []
chr_del = []
id_del = []
loc_lof = []
chr_lof = []
id_lof = []
pop_miss = []
pop_del = []
pop_lof = []

#%%
os.chdir('/home/prm/Desktop/Ptarmigan_analysis/vep_stuff/final/')
for file in glob.glob("Py_ind*_chrz.txt"):
    append_mutations(file, 'Py')
    append_geneID(file, 'Py')
    append_location(file, 'Py')
    save_geneID(file, 'Py')

for file in glob.glob("Alps_ind*_chrz.txt"):
    append_mutations(file, 'Alps')
    append_geneID(file, 'Alps')
    append_location(file, 'Alps')
    save_geneID(file, 'Alps')
        
for file in glob.glob("Ice_ind*_chrz.txt"):
    append_mutations(file, 'Ice')
    append_geneID(file, 'Ice')
    append_location(file,'Ice')
    save_geneID(file, 'Ice')
        
for file in glob.glob("Jam_ind*_chrz.txt"):
    append_mutations(file, 'Jam')
    append_geneID(file, 'Jam')
    append_location(file, 'Jam')
    save_geneID(file, 'Jam')

for file in glob.glob("Gr_ind*_chrz.txt"):
    append_mutations(file, 'GrE')
    append_geneID(file, 'GrE')
    append_location(file, 'GrE')
    save_geneID(file, 'GrE')
       
for file in glob.glob("GrW_ind*_chrz.txt"):
    append_mutations(file, 'GrW')
    append_geneID(file, 'GrW')
    append_location(file, 'GrW')
    save_geneID(file, 'GrW')

for file in glob.glob("Sva_ind*_chrz.txt"):
    append_mutations(file, 'Sva')
    append_geneID(file, 'Sva')
    append_location(file, 'Sva')
    save_geneID(file, 'Sva')

for file in glob.glob("NL_ind*_chrz.txt"):
    append_mutations(file, 'NL')
    append_geneID(file, 'NL')
    append_location(file, 'NL')
    save_geneID(file, 'NL')

for file in glob.glob("Scot_ind*_chrz.txt"):
    append_mutations(file, 'Scot')
    append_geneID(file, 'Scot')
    append_location(file, 'Scot')
    save_geneID(file, 'Scot')

for file in glob.glob("Lago_ind*_chrz.txt"):
    append_mutations(file, 'Lago')
    append_geneID(file, 'Lago')
    append_location(file, 'Lago')
    save_geneID(file, 'Lago')


# %%
vep_stats_columns = ['ind', 'pop', 'syn', 'missense', 'del', 'stop_gained', 'miss/syn', 'del/syn', 'LOF/syn']
vep_stats = pd.DataFrame(zip(ind, pop, syn_mut, miss_mut, delet_mut, lof_mut, miss_syn, delet_syn, lof_syn), columns = vep_stats_columns)
vep_stats.to_csv('vep_summary_stats_chrZ.txt', sep = '\t', header = True, index = False)

#%%
locations_columns = ['pop', 'ind', 'chr', 'location']
missense_mutations_locations = pd.DataFrame(zip(pop_miss, id_miss, chr_miss, loc_miss), columns = locations_columns)
missense_mutations_locations.to_csv('missense_locations_chrZ.txt', sep = '\t', header = True, index = False)
del_mutations_location = pd.DataFrame(zip(pop_del, id_del, chr_del, loc_del), columns = locations_columns)
del_mutations_location.to_csv('del_locations_chrZ.txt', sep = '\t', header = True, index = False)
lof_mutations_location = pd.DataFrame(zip(pop_lof, id_lof, chr_lof, loc_lof), columns = locations_columns)
lof_mutations_location.to_csv('lof_locations_chrZ.txt', sep = '\t', header = True, index = False)
