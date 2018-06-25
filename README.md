# mc_datasets

Bioinformatic analyses of amplicon and shotgun datasets used for investigating
MGS bias.

## Pre-requisites

### Dependencies

#### R packages

* `tidyverse`
* `dotenv` <https://github.com/gaborcsardi/dotenv>
* `rentrez` is needed only if re-generating the sample metadata files
* `dada2`
* `phyloseq`
* `metaphlanr` if planning to post-process Metaphlan2 output to a read table

#### Software

* NCBI's SRA Toolkit for downloading sequence data
* Metaphlan2 for microbiome profiling of shotgun data
* (optional) Aspera's `ascp` program for faster downloading of sequence data

### Environment file

A file called `.env` must be set up with various directory paths, which can be
parsed by the R package `dotenv`. An example is in `.env_example`.

## Datasets

### Amplicon

#### brooks2015

Brooks JP, Edwards DJ, Harwich MD, Rivera MC, Fettweis JM, Serrano MG, Reris
RA, Sheth NU, Huang B, Girerd P, Strauss JF, Jefferson KK, Buck GA. 2015. The
truth about metagenomics: quantifying and counteracting bias in 16S rRNA
studies. BMC Microbiol 15:66.

#### sinha2015

Sinha R, Abu-Ali G, Vogtmann E, Fodor AA, Ren B, Amir A, Schwager E, Crabtree
J, Ma S, Abnet CC, Knight R, White O, Huttenhower C. 2017. Assessment of
variation in microbial community amplicon sequencing by the Microbiome Quality
Control (MBQC) project consortium. Nat Biotechnol 35:1077–1086.

### Shotgun

#### costea2017

Costea PI, Zeller G, Sunagawa S, Pelletier E, Alberti A, Levenez F, Tramontano
M, Driessen M, Hercog R, Jung F-E, Kultima JR, Hayward MR, Coelho LP,
Allen-Vercoe E, Bertrand L, Blaut M, Brown JRM, Carton T, Cools-Portier S,
Daigneault M, Derrien M, Druesne A, de Vos WM, Finlay BB, Flint HJ, Guarner F,
Hattori M, Heilig H, Luna RA, van Hylckama Vlieg J, Junick J, Klymiuk I,
Langella P, Le Chatelier E, Mai V, Manichanh C, Martin JC, Mery C, Morita H,
O’Toole PW, Orvain C, Patil KR, Penders J, Persson S, Pons N, Popova M, Salonen
A, Saulnier D, Scott KP, Singh B, Slezak K, Veiga P, Versalovic J, Zhao L,
Zoetendal EG, Ehrlich SD, Dore J, Bork P. 2017. Towards standards for human
fecal sample processing in metagenomic studies. Nat Biotechnol 35:1069–1076.
