# FuncAnnot

## Descriptions

This pipeline is set to to non-model plant GO/KO function annotations using R. 

    I have to say it's a nightmare to blastp/blastx against NCBI NR protein database (gzipped 80G?) or find the domain using interproscan. Really a pain of patience.

Now it's good to try [Eggnog-mapper](http://eggnog-mapper.embl.de/) and swissprot, and refseq is another choice.

It can collect data from Eggnog-mapper of custom BLAST. But you may need to transform the seqID in database to the required format as described below.

- [x] Support Eggnog output '.annotations' files
- [x] Support custom GO/KO files

## Requirements

  **Perl modules**:

    JSON, Getopt::Long, [FuhaoPerl5Lib](https://github.com/lufuhao/FuhaoPerl5Lib)

  **R packages**: 

    optparse, dplyr, stringr, purrr, tidyverse

    AnnotationForge, ClusterProfiler, DBI, GO.db, RSQLite, enrichplot

  **Python**

    gzip, getopt, sys



## Options

|  Options     |  Parameters | Descriptions                                         |
| :--------    | :--------:  | :--------                                            |
| --help/-h    | -           | Print this help/usage                                |
| --prefix/-p  | STR         | Output prefix, default: MyOut                        |
|              |             |                                                      |
| #  For eggnog files                                                               |
| --input/-i   | Files       | Eggnog output files ending with ".annotations"       |
|              |             | Could be a list separated by comma or                |
|              |             | specify multiple times or both                       |
| --ignoreegg  | INT         | Ignore first INT comment lines, default:4            |
|              |             | Will inore lines starting with '# ' or '#query_name' |
| --rmeggtsnum |             | Trim transcript number ending with .\\d+$            |
|              |             |                                                      |
|#  For custom GOs                                                                  |
| --gos/-g     | Files       | 2-column data:                                       |
|              |             | column1: gene/transcript ID                          |
|              |             | column2: GOs sep by comma: GO:000000,GO000001        |
| --ignorego   | INT         | Ignore first INT comment lines, default:0            |
| --rmgotsnum  |             | Trim transcript number ending with .\\d+$            |
|              |             |                                                      |
| For custom KOs                                                                    |
| --kos/-k     | Files       | 2-column data:                                       |
|              |             |    column1: gene/transcript ID                       |
|              |             |    column2: GOs sep by comma. example: K00001,K00002 |
| --ignoreko   | INT         | Ignore first INT comment lines, default:0            |
| --rmkotsnum  |             | Trim transcript number ending with .\\d+$            |
| --keggjson   |             | KEGG Orthology file, downloaded from [KEGG](https://www.genome.jp/kegg-bin/show_brite?ko00001.keg) |


## Install

```
git clone https://github.com/lufuhao/FuncAnnot.git
export PATH=/your/path/to/FuncAnnot/bin:$PATH
```

## Steps:

### 1. prepare the data

#### 1.1. For eggnog

* Go to [Eggnog](http://eggnog-mapper.embl.de/) websites and upload your protein/DNA sequences for analysis
* Do fill the Email addess
* Change the parameters if necessary
* Submit
* An Email will be received, open the link and click 'start' to run the job
* Another Email will be received once the job is done; Click the link in Email, and download the results: .annotations .seed_orthologs (Optional) and info (Optional)

#### 1.2 For custom BLAST/Diamond/GhostX/RapSearch

* Run Diamond against Swissprot DB, get tab-delimited out
* Download IDmapping DB from [PIR](ftp://ftp.pir.georgetown.edu/databases/idmapping/idmapping.tb.gz)
```
wget ftp://ftp.pir.georgetown.edu/databases/idmapping/idmapping.tb.gz
```
* Transform swissID to GO using script/uniport2go.py
```
uniport2go.py -f res.blast.tab -m idmapping.tb.gz -o blast2go.tab
```

### 2. Data transform for AnnotationForge

```
cd $runningdir

funcannot.1.data.pl \
 -i 0.eggnog/HC001.annotations \
 -i 0.eggnog/HC002.annotations \
 -i 0.eggnog/stringtie.new.00000001.annotations \
 -i 0.eggnog/stringtie.new.00000002.annotations \
 -i 0.eggnog/stringtie.new.00000003.annotations \
 -i 0.eggnog/stringtie.new.00000004.annotations \
 -i 0.eggnog/LC001.annotations \
 -i 0.eggnog/LC002.annotations \
 -i 0.eggnog/LC003.annotations \
 --ignoreegg 4 --rmeggtsnum \
 --gos 0.diamond/HC01.gene2go \
 --gos 0.diamond/gene_go/HC02.gene2go \
 --gos 0.diamond/gene_go/LC01.gene2go \
 --gos 0.diamond/gene_go/LC02.gene2go \
 --gos 0.diamond/gene_go/LC03.gene2go \
 --gos 0.diamond/gene_go/stringtie.new.gene2go \
 --ignorego 0 --rmgotsnum -p MyDB
 --keggjson ko00001.json

head -n 3 *.tab
```


### 3. make custom OrgDB using R AnnotationForge

Do some research before running
* Go [NCBI](https://www.ncbi.nlm.nih.gov/taxonomy) and find your taxID, forexample, "4565" for Triticum aestivum
* Find the genus name: "Triticum"
* find and species name: "aestivum"
* define a version name: "\d+.\d+", otherwise version error would be reported
* Use the same prefix as Step2, for instance: **MyDB**
* It takes time to load all the data to DB depending on your data size, so pay patience ...

```bash
cd $runningdir

funcannot.2.makeOrgDB.rscript --prefix MyDB --genus Triticum --species aestivum --taxid 4565 --vers 2021.0128
```

```Rscript
install.packages("org.Taestivum.eg.db", repos = NULL, type = "source")
library(org.Taestivum.eg.db)
columns(org.Taestivum.eg.db)
keys(org.Taestivum.eg.db)
select(org.Taestivum.eg.db, keys="xxx", columns=c("GO, "Pathway"))
```

### 4. Given a gene list, do the enrichment of GO/KO using clusterProfiler


```
funcannot.3.annot.rscript --file "DEG.list" --prefix MyOut --lib "org.Taestivum.eg.db" --width 17 --height 17 --resolution 600 --pointsize 10 --family "ArialMT" 
```


## Revision history

v20211004 support emapper v1



## Author:

    Fu-Hao Lu

    Professor, PhD

    State Key Labortory of Crop Stress Adaptation and Improvement

    College of Life Science

    Jinming Campus, Henan University

    Kaifeng 475004, P.R.China

    E-mail: LUFUHAO@HENU.EDU.CN
