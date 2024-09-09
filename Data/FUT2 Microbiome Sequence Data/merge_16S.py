## Script to parse and combine the QIIME2 and VALENCIA outputs
# Used to combine with FUT2 clinical data and create the multiple sheets in FUT2Dataset.xlsx
import pandas as pds
import numpy as np

# Read Species level tables for each primerset
v1mixedv2 = pds.read_csv('./Species_table_V1MixedV2.tsv', delimiter='\t', skiprows=1)
v1v3 = pds.read_csv('./Species_table_V1V3.tsv', delimiter='\t', skiprows=1)
v1v2 = pds.read_csv('./Species_table_V1V2.tsv', delimiter='\t', skiprows=1)

v1v3_v1v2merge = v1v3.merge(v1v2, on='#OTU ID', how='outer')

allMerge = v1mixedv2.merge(v1v3_v1v2merge, on='#OTU ID', how='outer')

# Read Genera level tables for each primerset
v1mixedv2_Genera = pds.read_csv('./Genera_table_V1MixedV2.tsv', delimiter='\t', skiprows=1)
v1v3_Genera = pds.read_csv('./Genera_table_V1V3.tsv', delimiter='\t', skiprows=1)
v1v2_Genera = pds.read_csv('./Genera_table_V1V2.tsv', delimiter='\t', skiprows=1)

v1v3_v1v2merge_Genera = v1v3_Genera.merge(v1v2_Genera, on='#OTU ID', how='outer')
allMerge_Genera = v1mixedv2_Genera.merge(v1v3_v1v2merge_Genera, on='#OTU ID', how='outer')

# Read the FUT2 study data
FUT2Meta = pds.read_excel('../FUT2-Dataset.xlsx', sheet_name=0)

#
# Parse and merge the species matrix
# 
speciesMat = allMerge.iloc[:, 1:].T
speciesMat.columns = allMerge['#OTU ID'].values
speciesMat[pds.isna(speciesMat)] = 0

speciesMat = speciesMat.reset_index()
speciesMat = speciesMat.rename(columns={'index':'seqrun_samplenum'})
speciesMat['seqrun_samplenum'] = speciesMat['seqrun_samplenum'].str.replace('-MS28F-388R', '')
speciesMat['seqrun_samplenum'] = speciesMat['seqrun_samplenum'].str.replace('-28FCombo', '')
speciesMat['seqrun_samplenum'] = speciesMat['seqrun_samplenum'].str.replace('-MS28F', '')

# Calculate relative abundance to remove rare taxa (keep only taxa with counts > 0.05% relative abundance in  >= 2 samples). 
speciesMatRealAbundance = speciesMat.iloc[:, 1:].divide(speciesMat.iloc[:, 1:].sum(axis=1), axis=0)

speciesToRemove = speciesMatRealAbundance.columns[((speciesMatRealAbundance > 0.05).sum(axis=0) < 2)]

speciesMat_Filtered = speciesMat.drop(columns=speciesToRemove)

# curate the taxa name
def getLowestTaxonLevel_Species(taxaString):
    splitNameStrings = taxaString.split(';')
    currentLevel = -1
    currentLevelString = splitNameStrings[currentLevel]
    while currentLevelString == '__':
        currentLevel -= 1
        currentLevelString = splitNameStrings[currentLevel]
    if currentLevelString in ['s__uncultured_bacterium', 's__unidentified', 's__uncultured_Coriobacteriaceae']:
    	currentLevelString = splitNameStrings[currentLevel - 1] + currentLevelString[2:]
    return currentLevelString
    
speciesNamesParsed = [getLowestTaxonLevel_Species(x) for x in speciesMat_Filtered.columns]

# Remove the s__ from species
speciesNamesParsed = [x.replace('s__', '') for x in speciesNamesParsed]
# Remove the g__ from genera
speciesNamesParsed = [x.replace('g__', '') for x in speciesNamesParsed]
# Remove the g__ from genera
speciesNamesParsed = [x.replace('d__', 'Domain ') for x in speciesNamesParsed]
speciesNamesParsed = [x.replace('f__', 'Family ') for x in speciesNamesParsed]
speciesNamesParsed = [x.replace('c__', 'Class ') for x in speciesNamesParsed]

speciesNamesParsed = [x.replace('_', ' ') if x!= 'seqrun_samplenum' else x for x in speciesNamesParsed]

speciesMat_Filtered.columns = speciesNamesParsed

speciesDataset = FUT2Meta.merge(speciesMat_Filtered, on='seqrun_samplenum')

speciesDataset.to_csv('./FUT2_SpeciesMatrix.csv', index=False)
#
# Parse and merge the genera dataset
#

generaMat = allMerge_Genera.iloc[:, 1:].T
generaMat.columns = allMerge_Genera['#OTU ID'].values
generaMat[pds.isna(generaMat)] = 0

generaMat = generaMat.reset_index()
generaMat = generaMat.rename(columns={'index':'seqrun_samplenum'})
generaMat['seqrun_samplenum'] = generaMat['seqrun_samplenum'].str.replace('-MS28F-388R', '')
generaMat['seqrun_samplenum'] = generaMat['seqrun_samplenum'].str.replace('-28FCombo', '')
generaMat['seqrun_samplenum'] = generaMat['seqrun_samplenum'].str.replace('-MS28F', '')

# Calculate relative abundance to remove rare taxa (keep only taxa with counts > 0.05% relative abundance in  >= 2 samples). 
generaMatRealAbundance = generaMat.iloc[:, 1:].divide(generaMat.iloc[:, 1:].sum(axis=1), axis=0)

generaToRemove = generaMatRealAbundance.columns[((generaMatRealAbundance > 0.05).sum(axis=0) < 2)]

generaMat_Filtered = generaMat.drop(columns=generaToRemove)

# Recalculate the relative abundance with the rare taxa filtered
generaMatRealAbundance = generaMat_Filtered.iloc[:, 1:].divide(generaMat_Filtered.iloc[:, 1:].sum(axis=1), axis=0)

generaMatRealAbundance = pds.concat([generaMat_Filtered.iloc[:, 0], generaMatRealAbundance], axis=1)


# curate the taxa name
def getLowestTaxonLevel_Genera(taxaString):
    splitNameStrings = taxaString.split(';')
    currentLevel = -1
    currentLevelString = splitNameStrings[currentLevel]
    while currentLevelString == '__':
        currentLevel -= 1
        currentLevelString = splitNameStrings[currentLevel]
    return currentLevelString
    
generaNamesParsed = [getLowestTaxonLevel_Genera(x) for x in generaMat_Filtered.columns]

# Remove the g__ from genera
generaNamesParsed = [x.replace('g__', '') for x in generaNamesParsed]
# Remove the g__ from genera
generaNamesParsed = [x.replace('d__', 'Domain ') for x in generaNamesParsed]
generaNamesParsed = [x.replace('f__', 'Family ') for x in generaNamesParsed]
generaNamesParsed = [x.replace('c__', 'Class ') for x in generaNamesParsed]

generaMat_Filtered.columns = generaNamesParsed

generaMatRealAbundance.columns = generaNamesParsed
 
generaDataset = FUT2Meta.merge(generaMatRealAbundance, on='seqrun_samplenum')
generaDataset.to_csv('./FUT2_GeneraMatrix.csv', index=False)


# Combine the ASV tables
v1v2_asv = pds.read_csv('./table_asv_V1V2.tsv', delimiter='\t', skiprows=1)
v1v2_asv_taxonomy = pds.read_csv('./taxonomy_V1V2.tsv', delimiter='\t')

v1mixedv2_asv = pds.read_csv('./table_asv_V1MixedV2.tsv', delimiter='\t', skiprows=1)
v1mixedv2_asv_taxonomy = pds.read_csv('./taxonomy_V1MixedV2.tsv', delimiter='\t')

v1v3_asv = pds.read_csv('./table_asv_V1V3.tsv', delimiter='\t', skiprows=1)
v1v3_asv_taxonomy = pds.read_csv('./taxonomy_V1V3.tsv', delimiter='\t')

v1v3_v1v2merge_asv = v1v3_asv.merge(v1v2_asv, on='#OTU ID', how='outer')
allMerge_asv = v1mixedv2_asv.merge(v1v3_v1v2merge_asv, on='#OTU ID', how='outer')
allMerge_asv = allMerge_asv.fillna(0)


asvMat = allMerge_asv.iloc[:, 1:].T
asvMat.columns = allMerge_asv['#OTU ID'].values

asvMat = asvMat.reset_index()
asvMat = asvMat.rename(columns={'index':'seqrun_samplenum'})
asvMat['seqrun_samplenum'] = asvMat['seqrun_samplenum'].str.replace('-MS28F-388R', '')
asvMat['seqrun_samplenum'] = asvMat['seqrun_samplenum'].str.replace('-28FCombo', '')
asvMat['seqrun_samplenum'] = asvMat['seqrun_samplenum'].str.replace('-MS28F', '')

asvDataset = FUT2Meta.merge(asvMat, on='seqrun_samplenum')

asvDataset = pds.concat([asvDataset['SampleID'], asvDataset.iloc[:, 35:]], axis=1)
asvDataset = asvDataset.T

asvDataset_SampleNamesColumns = asvDataset.iloc[0, :]

asvDataset = asvDataset.iloc[1:, :]
asvDataset.columns = asvDataset_SampleNamesColumns
asvDataset.rename(index={'SampleID':'ASVID'}, inplace=True)

asvDataset.to_csv('./FUT2_asvMatrix.csv', index=True)

# Combine the ASV taxonomy
v1v3_v1v2merge_asv_taxonomy = v1v3_asv_taxonomy.merge(v1v2_asv_taxonomy, on='Feature ID', how='outer')
v1v3_v1v2merge_asv_taxonomy = v1v3_v1v2merge_asv_taxonomy.drop(columns=['Confidence_x', 'Confidence_y'])
v1v3_v1v2merge_asv_taxonomy.loc[pds.isna(v1v3_v1v2merge_asv_taxonomy['Taxon_x']), 'Taxon_x'] = v1v3_v1v2merge_asv_taxonomy.loc[pds.isna(v1v3_v1v2merge_asv_taxonomy['Taxon_x']), 'Taxon_y']
v1v3_v1v2merge_asv_taxonomy = v1v3_v1v2merge_asv_taxonomy.drop(columns=['Taxon_y'])
v1v3_v1v2merge_asv_taxonomy.rename(columns={'Taxon_x':'Taxon'}, inplace=True)

allMerge_asv_taxonomy = v1mixedv2_asv_taxonomy.merge(v1v3_v1v2merge_asv_taxonomy, on='Feature ID', how='outer')
allMerge_asv_taxonomy = allMerge_asv_taxonomy.drop(columns=['Confidence'])
allMerge_asv_taxonomy.loc[pds.isna(allMerge_asv_taxonomy['Taxon_x']), 'Taxon_x'] = allMerge_asv_taxonomy.loc[pds.isna(allMerge_asv_taxonomy['Taxon_x']), 'Taxon_y']
allMerge_asv_taxonomy = allMerge_asv_taxonomy.drop(columns=['Taxon_y'])
allMerge_asv_taxonomy.rename(columns={'Taxon_x':'Taxon'}, inplace=True)

# allMerge_asv_taxonomy.to_csv('./FUT2_asvTaxonomy.csv', index=False)

# Valencia parse the results
v1mixedv2_Valencia = pds.read_csv('./VALENCIA_V1MixedV2.csv')
v1v2_Valencia = pds.read_csv('./VALENCIA_V1V2.csv')
v1v3_Valencia = pds.read_csv('./VALENCIA_V1V3.csv')

v1mixedv2_Valencia.rename(columns={'sampleID':'seqrun_samplenum'}, inplace=True)
v1v2_Valencia.rename(columns={'sampleID':'seqrun_samplenum'}, inplace=True)
v1v3_Valencia.rename(columns={'sampleID':'seqrun_samplenum'}, inplace=True)

v1mixedv2_Valencia = v1mixedv2_Valencia.loc[:, ['seqrun_samplenum', 'subCST', 'score', 'CST']]
v1v2_Valencia = v1v2_Valencia.loc[:, ['seqrun_samplenum', 'subCST', 'score', 'CST']]
v1v3_Valencia = v1v3_Valencia.loc[:, ['seqrun_samplenum', 'subCST', 'score', 'CST']]
valenciaAll = pds.concat([v1mixedv2_Valencia, v1v2_Valencia, v1v3_Valencia], axis=0)

valenciaAll['seqrun_samplenum'] = valenciaAll['seqrun_samplenum'].str.replace('-MS28F-388R', '')
valenciaAll['seqrun_samplenum'] = valenciaAll['seqrun_samplenum'].str.replace('-28FCombo', '')
valenciaAll['seqrun_samplenum'] = valenciaAll['seqrun_samplenum'].str.replace('-MS28F', '')
valenciaAll.rename(columns={'score':'VALENCIA_score'}, inplace=True)

# drop columns if they exist already
FUT2Meta.drop(columns=['subCST', 'VALENCIA_score', 'CST', 'LactobacillusDominant_CST', 'LactobacillusDominant'], inplace=True)

# lactDom = generaDataset.loc[:, 34:].apply(np.argmax, axis=1)
lactDom = generaDataset.loc[:, 'Lactobacillus'] > 0.9
lactDom[lactDom == True] = 'L. Dominant'
lactDom[lactDom == False] = 'L. Depleted'

LactoDomGeneraRelAbundance = pds.concat([generaDataset.loc[:, ['seqrun_samplenum']], lactDom], axis=1)

LactoDomGeneraRelAbundance.rename(columns={'Lactobacillus':'LactobacillusDominant'}, inplace=True)

valenciaAll = valenciaAll.merge(LactoDomGeneraRelAbundance, on='seqrun_samplenum')

# Not required as Lactobacillus dominant CST based criteria can also be set up in the R scripts
LactobacillusDominant_CST =  valenciaAll['CST'].str.contains('IV')
LactobacillusDominant_CST[LactobacillusDominant_CST == False] = 'L. Dominant'
LactobacillusDominant_CST[LactobacillusDominant_CST == True] = 'L. Depleted'

valenciaAll.insert(valenciaAll.shape[1], 'LactobacillusDominant_CST' , LactobacillusDominant_CST)

FUT2Meta = FUT2Meta.merge(valenciaAll, on='seqrun_samplenum')

FUT2Meta.to_csv('FUT2 Metadata with CST.csv', index=False)
