# FuncAnnot

## Descriptions

This pipeline is set to to non-model plant GO/KO function annotations using R. 

    I have to say it's a nightmare to blastp/blastx against NCBI NR protein database (gzipped 80G?) or find the domain using interproscan. Really a pain of patience.

Now it's good to try [Eggnog-mapper](http://eggnog-mapper.embl.de/) and swissprot, and refseq is another choice.

It can collect data from Eggnog-mapper of custom BLAST. But you may need to transform the seqID in database to the required format as described below.

[x] Support Eggnog output '.annotations' files

[x] Support custom GO/KO files


## Options

| Options   |  Parameters | Descriptions |
| --------  | --------  | --------  |
| --help/-h | - | Print this help/usage |
| --prefix/-p | STR | Output prefix, default: MyOut |
| #  For eggnog files
| --input|-i | Files | Eggnog output files ending with ".annotations" |

>                           Could be a list separated by comma or

>                           specify multiple times or both

>    --ignoreegg   <INT>    Ignore first INT comment lines, default:4

>                           Will inore lines starting with '# ' or '#query_name'

>    --rmeggtsnum           Trim transcript number ending with .\\d+$

>

>#  For custom GOs

>    --gos|-g      <Files> 2-column data: 

>                              column1: gene/transcript ID

>                              column2: GOs sep by comma. example: GO:000000,GO000001

>    --ignorego    <INT>   Ignore first INT comment lines, default:0

>    --rmgotsnum           Trim transcript number ending with .\\d+$

>

>#  For custom KOs

>    --kos|-k      <Files> 2-column data: 

>                              column1: gene/transcript ID

>                              column2: GOs sep by comma. example: K00001,K00002

>    --ignoreko    <INT>   Ignore first INT comment lines, default:0

>    --rmkotsnum           Trim transcript number ending with .\\d+$

>    --keggjson            KEGG Orthology file, downloaded from [KEGG](https://www.genome.jp/kegg-bin/show_brite?ko00001.keg)







## Author:

    Fu-Hao Lu

    Professor, PhD

    State Key Labortory of Crop Stress Adaptation and Improvement

    College of Life Science

    Jinming Campus, Henan University

    Kaifeng 475004, P.R.China

    E-mail: LUFUHAO@HENU.EDU.CN
