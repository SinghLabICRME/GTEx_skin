## Project: Analysis of sun-exposed skin transcriptome between slow and fast death cases
#### Author: "Ahmed S Abouhashem"
#### Date: "22/03/2022"

### Download data related files from GTEx Portal using command line
Download gene TPM values
```
wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
```
Download sample attributes
```
wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
```
Download description of sample attributes (RNA quality, tissue source, …etc)
```
wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx
```
Download subject phenotypes
```
wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx
```
Download Description of subject phenotypes: gender, Age and hardy scale (death type: slow, fast violent, fast unexpected, intermediate, ventilator or unknown)
```
wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
```

### Filter sun-exposed skin samples and their metadata using Python
Read the gene TPM values
```
import pandas as pd
all_samples = pd.read_csv('/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct', sep="\t",skiprows=2)
all_samples = all_samples.set_index('Name')
len(list(all_samples.columns)) #17384 (total number of samples)
```
Drop the description column to keep gene IDs and TPM values only
```
all_samples = all_samples.drop(columns = 'Description')
```
Read the metadata
```
metadata = pd.read_csv('/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', sep="\t")
```
Select only sun-exposed skin samples with available RNA-seq data
```
skin_samples_metadata = metadata.loc[(metadata['SMTSD'] == 'Skin - Sun Exposed (Lower leg)') & (metadata['SMAFRZE'] == 'RNASEQ')]
skin_samples_IDs = skin_samples_metadata[skin_samples_metadata.columns[0]]
```
Subset sun exposed skin tpm values from all counts matrix
```
skin_samples_tpm = pd.DataFrame(all_samples, columns=skin_samples_IDs)
```
Extract skin samples TPM matrix and metadata to as a csv file
```
skin_samples_tpm.to_csv(r'/Skin_samples_tpm.csv')
skin_samples_metadata = skin_samples_metadata.set_index('SAMPID')
skin_samples_metadata.to_csv(r'/Skin_samples_metadata.csv')
```
### Add phenotypes data (AGE, SEX and death type) to the original metadata file
First, convert both samples metadata and phenotypes data to excel files (*.xlsx). Next, create a dictionary of individuals phenotypes data.
```
import openpyxl
filepath = '/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.xlsx'
wb = openpyxl.load_workbook(filepath, data_only=True)
sheet = wb.active
c = sheet.cell

phenotypes = {}
for i in range(2,982):
    key = c(row=i,column=1).value
    value = []
    for j in range(2,5):
        value.append(c(row=i,column=j).value)
    phenotypes[key]=value
```
Next, Add 5 columns to the original metadata file (Individual ID, AGE, SEX, DTHHRDY (hardy.scale) and death_type)
```
filepath = '/mnt/c87fe4a6-44ed-460c-b6a1-d4ebbce31384/Temp/Kanhaiya meeting/Post mortem project/2022/Final/Skin_samples_metadata.xlsx'
wb = openpyxl.load_workbook(filepath, data_only=True)
sheet = wb.active
c = sheet.cell

#Write the new column headers
c(row=1,column=64).value = 'individual.id'
c(row=1,column=65).value = 'SEX'
c(row=1,column=66).value = 'AGE'
c(row=1,column=67).value = 'DTHHRDY'
c(row=1,column=68).value = 'death_type'

#Fetch sample IDs from each sample name (example: this sample ID "GTEX-1PIIG-0226-SM-EWRN6" includes "GTEX-1PIIG" as individual ID). This part is to match with the phenotypes IDs in the dictionary "phenotypes" created before.
for i in range(2,703):
    current_sample = c(row=i,column=1).value
    dash = 0
    id = ''
    for letter in current_sample:
        if letter == '-':
            dash += 1
            if dash == 2:
                break
            else:
                id += letter
        else:
            id += letter
    c(row=i, column=64).value = id
    c(row=i,column=65).value = phenotypes[id][0]
    c(row=i, column=66).value = phenotypes[id][1]
    c(row=i, column=67).value = phenotypes[id][2]

    #Classify samples based on the provided hardy.scale
    if c(row=i, column=67).value == None:
        c(row=i, column=68).value = 'Unknown'
    if c(row=i, column=67).value == 0:
        c(row=i, column=68).value = 'Ventilator'
    if c(row=i, column=67).value == 1:
        c(row=i, column=68).value = 'Fast Violent'
    if c(row=i, column=67).value == 2:
        c(row=i, column=68).value = 'Fast Unexpected'
    if c(row=i, column=67).value == 3:
        c(row=i, column=68).value = 'Intermediate'
    if c(row=i, column=67).value == 4:
        c(row=i, column=68).value = 'Slow'

#Overwrite the existing file and then convert it as a csv file again.
wb.save(filepath)
```
### Using R: Create a Seurat object including sun exposed skin samples
Import the counts matrix of sun-exposed skin samples and their metadata metadata
```
counts = read.csv('/Skin_samples_tpm.csv', sep = ',',row.names = 1)
metadata = read.csv('/combined_metadata.csv', sep = ',',row.names = 1)
```
subset '-' with '.' in metadata row names to match with column names
```
rownames(metadata) = gsub('-','.',rownames(metadata))
```
Create a seurat object
```
object = CreateSeuratObject(counts = counts, meta.data = metadata,project = 'GTEX-1RB15 GTEx Analysis V8 (dbGaP Accession phs000424.v8.p2)')
```
Data log normalization
```
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 1000000)
```
Data log normalization
```
Differential expression analysis
```
Idents(object)='death_type'
markers = FindMarkers(object, ident.1 = 'Slow',ident.2 = c('Fast Violent','Fast Unexpected'),min.pct = 0.2,logfc.threshold = 0)
write.csv(markers, "Supplementary file 2 (DEGs).csv")
```
### Next, perform the enrichment analysis using g:Profiler using web interface for upregulated and downregulated genes with adjusted p value < 0.01 and log2FC ±0.5 (https://biit.cs.ut.ee/gprofiler/gost).
