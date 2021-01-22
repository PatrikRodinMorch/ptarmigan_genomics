#!/usr/bin/python

import pandas as pd 
import argparse
import allel

parser = argparse.ArgumentParser(description = 'Script that returns annotations for ROHs and overlap between ROHs and vep mutation classes')
parser.add_argument('-pop', type = str, help = 'The population under consideration', required = True)

opts = parser.parse_args()

#%%
annotations = allel.gff3_to_dataframe("/home/prm/Desktop/Ptarmigan_analysis/ROH/Gallus_gallus.GRCg6a.99.chr.gff3", attributes = ['type', 'start', 'end', 'ID', 'biotype', 'description'])


def roh_annotation(pop, annotations):
    plink = pd.read_table('/home/prm/Desktop/Ptarmigan_analysis/ROH/'+pop+'/plink.hom', sep = '\s+')

    annotations_start = annotations.iloc[:, 3].tolist()
    annotations_end = annotations.iloc[:, 4].tolist()
    AID = annotations.iloc[:, 11].str.split(':', expand = True)
    annotations_ID = AID.iloc[:, 1].tolist()
    annotations_chrom = annotations.iloc[:, 0].tolist()
    annotations_type = annotations.iloc[:, 2].tolist()

    plink_start = plink['POS1'].tolist()
    plink_end = plink['POS2'].tolist()
    plink_ID = plink['IID'].tolist()
    plink_ROH = plink['KB'].tolist()
    plink_chrom = plink['CHR'].tolist()
    plink['new_chrom'] = 'Chr' + plink['CHR'].astype(str)
    plink_new_chr = plink['new_chrom'].tolist()



    ROH_anno = []

    for i in range(len(annotations_start)):
        for n in range(len(plink_start)):
            if (annotations_type[i] != 'chromosome') & (annotations_type[i] == 'gene'):
                if (annotations_chrom[i] == str(plink_chrom[n])) and (plink_start[n] <= annotations_start[i]) and (plink_end[n] >= annotations_end[i]):
                    annot = [annotations_chrom[i],annotations_start[i], annotations_end[i], annotations_type[i], annotations_ID[i], plink_chrom[n], plink_start[n], plink_end[n], plink_ID[n], plink_ROH[n]]
                    ROH_anno.append(annot)            

    columns = ['gff_chr', 'gff_start', 'gff_end', 'gff_type', 'gff_geneID', 'plink_chr', 'plink_start', 'plink_end', 'plink_id', 'plink_ROH_size']
    ROH_annotation = pd.DataFrame(ROH_anno, columns = columns)
   #ROH_annotation.to_csv('/home/prm/Desktop/Ptarmigan_analysis/ROH/'+pop+'/'+pop+'_ROH_annotation.txt', sep = '\t', header = True)
    ROH_annotation_uniq = ROH_annotation['gff_geneID'].unique()
    ROH_annotation_uniq.to_csv('/home/prm/Desktop/Ptarmigan_analysis/ROH/'+pop+'/'+pop+'_ROH_annotation_uniq.txt', sep = '\t', header = False)

    return plink_start, plink_end, plink_ID, plink_ROH, plink_new_chr, ROH_annotation


# %%

def missense_roh_overlap(pop):
    plink_start, plink_end, plink_ID, plink_ROH, plink_new_chr, ROH_annotation = roh_annotation(pop, annotations)
    
    mut = pd.read_table('/home/prm/Desktop/Ptarmigan_analysis/vep_stuff/final/missense_locations.txt', sep = '\t')

    mut_ind = mut['ind'].tolist()
    mut_pop = mut['pop'].tolist()
    mut_chr = mut['chr'].tolist()
    mut_loc = mut['location'].tolist()
    anno_ind = ROH_annotation['plink_id'].tolist()
    anno_chr = ROH_annotation['plink_chr'].tolist()
    anno_start = ROH_annotation['plink_start'].tolist()
    anno_end = ROH_annotation['plink_end'].tolist()
    anno_gene = ROH_annotation['gff_geneID'].tolist()
    ROH_annotation['new_chrom'] = 'Chr' + ROH_annotation['plink_chr'].astype(str)
    anno_new_chr = ROH_annotation['new_chrom'].tolist()
    

    mut_roh_overlap = []
    mut_in_ROH = 0
    total_number_mut = len(mut[mut['pop'] == pop])
    total_number_roh = len(plink_start)
    m_r_overlap = []
    m_r_overlap_geneID = []

    for i in range(len(plink_start)):
        for n in range(len(mut_loc)):
            if (str(mut_chr[n]) == str(plink_new_chr[i])) and (mut_ind[n] == plink_ID[i]) and (mut_loc[n] >= plink_start[i]) and (mut_loc[n] <= plink_end[i]):
                overlap = [mut_ind[n], mut_pop[n], mut_chr[n], mut_loc[n], plink_ID[i], plink_new_chr[i], plink_start[i], plink_end[i], plink_ROH[i]]
                m_r_overlap.append(overlap)
                mut_in_ROH += 1
            
    for i in range(len(mut_loc)):
        for n in range(len(anno_start)):
            if (str(mut_chr[i]) == str(anno_new_chr[n])) and (mut_ind[i] == anno_ind[n]) and (mut_loc[i] >= anno_start[n]) and (mut_loc[i] <= anno_end[n]):
                overlap_ID = [mut_ind[i], mut_pop[i], mut_chr[i], mut_loc[i], anno_ind[n], anno_new_chr[n], anno_start[n], anno_end[n], anno_gene[n]]
                m_r_overlap_geneID.append(overlap_ID)
      
    columns_overlap = ['mut_ind', 'mut_pop', 'mut_chr', 'mut_loc', 'plink_ID', 'plink_chrom', 'plink_start', 'plink_end', 'plink_roh']
    mut_roh_overlap = pd.DataFrame(m_r_overlap, columns = columns_overlap)
    #mut_roh_overlap.to_csv('/home/prm/Desktop/Ptarmigan_analysis/ROH/'+pop+'/'+pop+"_mut_missense_overlap.txt", sep = '\t', header = True)


    columns_overlap_ID = ['mut_ind', 'mut_pop', 'mut_chr', 'mut_loc', 'anno_ind', 'anno_chrom', 'anno_start', 'anno_end', 'anno_gene']
    mut_roh_overlap_geneID = pd.DataFrame(m_r_overlap_geneID, columns = columns_overlap_ID)
    #mut_roh_overlap_geneID.to_csv('/home/prm/Desktop/Ptarmigan_analysis/ROH/'+pop+'/'+pop+"_mut_missense_overlap_geneID.txt", sep = '\t', header = True)
    mut_roh_overlap_geneID_uniq = mut_roh_overlap_geneID['anno_gene'].unique()
    mut_roh_overlap_geneID_uniq.to_csv('/home/prm/Desktop/Ptarmigan_analysis/ROH/'+pop+'/'+pop+"_mut_missense_overlap_geneID_uniq.txt", sep = '\t', header = False)


    mut_stats = {'Population':pop,'Mutation inside ROH':mut_in_ROH, 'Mutations outside ROH':total_number_mut - mut_in_ROH, 'Percentage of mutations inside ROH': mut_in_ROH/total_number_mut , 'Total Number fo mutations':total_number_mut, 'Total number of ROH':total_number_roh}
    mut_stats_df = pd.DataFrame.from_dict(mut_stats, orient = 'index')
    #mut_stats_df.to_csv('/home/prm/Desktop/Ptarmigan_analysis/ROH/'+pop+'/'+pop+'_missense_stats.txt', sep = '\t', header = False)


def deleterious_roh_overlap(pop):
    plink_start, plink_end, plink_ID, plink_ROH, plink_new_chr, ROH_annotation = roh_annotation(pop, annotations)
    
    mut = pd.read_table('/home/prm/Desktop/Ptarmigan_analysis/vep_stuff/final/del_locations.txt', sep = '\t')


    mut_ind = mut['ind'].tolist()
    mut_pop = mut['pop'].tolist()
    mut_chr = mut['chr'].tolist()
    mut_loc = mut['location'].tolist()
    anno_ind = ROH_annotation['plink_id'].tolist()
    anno_chr = ROH_annotation['plink_chr'].tolist()
    anno_start = ROH_annotation['plink_start'].tolist()
    anno_end = ROH_annotation['plink_end'].tolist()
    anno_gene = ROH_annotation['gff_geneID'].tolist()
    ROH_annotation['new_chrom'] = 'Chr' + ROH_annotation['plink_chr'].astype(str)
    anno_new_chr = ROH_annotation['new_chrom'].tolist()

    mut_roh_overlap = []
    mut_in_ROH = 0
    total_number_mut = len(mut[mut['pop'] == pop])
    total_number_roh = len(plink_start)
    m_r_overlap = []
    m_r_overlap_geneID = []

    for i in range(len(plink_start)):
        for n in range(len(mut_loc)):
            if (str(mut_chr[n]) == str(plink_new_chr[i])) and (mut_ind[n] == plink_ID[i]) and (mut_loc[n] >= plink_start[i]) and (mut_loc[n] <= plink_end[i]):
                overlap = [mut_ind[n], mut_pop[n], mut_chr[n], mut_loc[n], plink_ID[i], plink_new_chr[i], plink_start[i], plink_end[i], plink_ROH[i]]
                m_r_overlap.append(overlap)

    for i in range(len(mut_loc)):
        for n in range(len(anno_start)):
            if (str(mut_chr[i]) == str(anno_new_chr[n])) and (mut_ind[i] == anno_ind[n]) and (mut_loc[i] >= anno_start[n]) and (mut_loc[i] <= anno_end[n]):
                overlap_ID = [mut_ind[i], mut_pop[i], mut_chr[i], mut_loc[i], anno_ind[n], anno_new_chr[n], anno_start[n], anno_end[n], anno_gene[n]]
                m_r_overlap_geneID.append(overlap_ID)

      
    columns_overlap = ['mut_ind', 'mut_pop', 'mut_chr', 'mut_loc', 'plink_ID', 'plink_chrom', 'plink_start', 'plink_end', 'plink_roh']
    mut_roh_overlap = pd.DataFrame(m_r_overlap, columns = columns_overlap)
   # mut_roh_overlap.to_csv('/home/prm/Desktop/Ptarmigan_analysis/ROH/'+pop+'/'+pop+"_mut_deleterious_overlap.txt", sep = '\t', header = True)

    columns_overlap_ID = ['mut_ind', 'mut_pop', 'mut_chr', 'mut_loc', 'anno_ind', 'anno_chrom', 'anno_start', 'anno_end', 'anno_gene']
    mut_roh_overlap_geneID = pd.DataFrame(m_r_overlap_geneID, columns = columns_overlap_ID)
    #mut_roh_overlap_geneID.to_csv('/home/prm/Desktop/Ptarmigan_analysis/ROH/'+pop+'/'+pop+"_mut_delet_overlap_geneID.txt", sep = '\t', header = True)
    mut_roh_overlap_geneID_uniq = mut_roh_overlap_geneID['anno_gene'].unique()
    mut_roh_overlap_geneID_uniq.to_csv('/home/prm/Desktop/Ptarmigan_analysis/ROH/'+pop+'/'+pop+"_mut_delet_overlap_geneID_uniq.txt", sep = '\t', header = False)



    mut_stats = {'Population':pop,'Mutation inside ROH':mut_in_ROH, 'Mutations outside ROH':total_number_mut - mut_in_ROH, 'Percentage of mutations inside ROH': mut_in_ROH/total_number_mut , 'Total Number fo mutations':total_number_mut, 'Total number of ROH':total_number_roh}
    mut_stats_df = pd.DataFrame.from_dict(mut_stats, orient = 'index')
    #mut_stats_df.to_csv('/home/prm/Desktop/Ptarmigan_analysis/ROH/'+pop+'/'+pop+'_deleterious_stats.txt', sep = '\t', header = False)

def lof_roh_overlap(pop):
    plink_start, plink_end, plink_ID, plink_ROH, plink_new_chr, ROH_annotation = roh_annotation(pop, annotations)
    
    mut = pd.read_table('/home/prm/Desktop/Ptarmigan_analysis/vep_stuff/final/lof_locations.txt', sep = '\t')


    mut_ind = mut['ind'].tolist()
    mut_pop = mut['pop'].tolist()
    mut_chr = mut['chr'].tolist()
    mut_loc = mut['location'].tolist()
    anno_ind = ROH_annotation['plink_id'].tolist()
    anno_chr = ROH_annotation['plink_chr'].tolist()
    anno_start = ROH_annotation['plink_start'].tolist()
    anno_end = ROH_annotation['plink_end'].tolist()
    anno_gene = ROH_annotation['gff_geneID'].tolist()
    ROH_annotation['new_chrom'] = 'Chr' + ROH_annotation['plink_chr'].astype(str)
    anno_new_chr = ROH_annotation['new_chrom'].tolist()

    mut_roh_overlap = []
    mut_in_ROH = 0
    total_number_mut = len(mut[mut['pop'] == pop])
    total_number_roh = len(plink_start)
    m_r_overlap = []
    m_r_overlap_geneID = []

    for i in range(len(plink_start)):
        for n in range(len(mut_loc)):
            if (str(mut_chr[n]) == str(plink_new_chr[i])) and (mut_ind[n] == plink_ID[i]) and (mut_loc[n] >= plink_start[i]) and (mut_loc[n] <= plink_end[i]):
                overlap = [mut_ind[n], mut_pop[n], mut_chr[n], mut_loc[n], plink_ID[i], plink_new_chr[i], plink_start[i], plink_end[i], plink_ROH[i]]
                m_r_overlap.append(overlap)
                mut_in_ROH += 1

    for i in range(len(mut_loc)):
            for n in range(len(anno_start)):
                if (str(mut_chr[i]) == str(anno_new_chr[n])) and (mut_ind[i] == anno_ind[n]) and (mut_loc[i] >= anno_start[n]) and (mut_loc[i] <= anno_end[n]):
                    overlap_ID = [mut_ind[i], mut_pop[i], mut_chr[i], mut_loc[i], anno_ind[n], anno_new_chr[n], anno_start[n], anno_end[n], anno_gene[n]]
                    m_r_overlap_geneID.append(overlap_ID)

      
    columns_overlap = ['mut_ind', 'mut_pop', 'mut:chr', 'mut_loc', 'plink_ID', 'plink_chrom', 'plink_start', 'plink_end', 'plink_roh']
    mut_roh_overlap = pd.DataFrame(m_r_overlap, columns = columns_overlap)
    #mut_roh_overlap.to_csv('/home/prm/Desktop/Ptarmigan_analysis/ROH/'+pop+'/'+pop+"_mut_lof_overlap.txt", sep = '\t', header = True)
    
    columns_overlap_ID = ['mut_ind', 'mut_pop', 'mut_chr', 'mut_loc', 'anno_ind', 'anno_chrom', 'anno_start', 'anno_end', 'anno_gene']
    mut_roh_overlap_geneID = pd.DataFrame(m_r_overlap_geneID, columns = columns_overlap_ID)
    #mut_roh_overlap_geneID.to_csv('/home/prm/Desktop/Ptarmigan_analysis/ROH/'+pop+'/'+pop+"_mut_lof_overlap_geneID.txt", sep = '\t', header = True)
    mut_roh_overlap_geneID_uniq = mut_roh_overlap_geneID['anno_gene'].unique()
    mut_roh_overlap_geneID_uniq.to_csv('/home/prm/Desktop/Ptarmigan_analysis/ROH/'+pop+'/'+pop+"_mut_lof_overlap_geneID_uniq.txt", sep = '\t', header = False)



    mut_stats = {'Population':pop,'Mutation inside ROH':mut_in_ROH, 'Mutations outside ROH':total_number_mut - mut_in_ROH, 'Percentage of mutations inside ROH': mut_in_ROH/total_number_mut , 'Total Number fo mutations':total_number_mut, 'Total number of ROH':total_number_roh}
    mut_stats_df = pd.DataFrame.from_dict(mut_stats, orient = 'index')
    #mut_stats_df.to_csv('/home/prm/Desktop/Ptarmigan_analysis/ROH/'+pop+'/'+pop+'_lof_stats.txt', sep = '\t', header = False)





roh_annotation(opts.pop, annotations)
missense_roh_overlap(opts.pop)
deleterious_roh_overlap(opts.pop)
lof_roh_overlap(opts.pop)
