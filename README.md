# hgdp
Package with data cleaned from the Human Genome Diversity Panel (HGDP) dataset. 

## Populations :

- Adygei in Russia
- Balochi in Pakistan
- Bandu in Kenya
- Basque in France
- Bedouin in Israel
- Bergamo in Italy
- Biaka in Central Africa Republic
- Bougainville in Bougainville
- Brahui in Pakistan
- Burusho in Pakistan
- Cambodian in Cambodia
- Colombian in Colombia
- Dai in China
- Daur in China
- Druze in Israel
- French in France
- Han in China
- Hazara in Pakistan
- Hezhen in China
- Japanese in Japan
- Kalash in Pakistan
- Karitiana in Brazil
- Lahu in China
- Makrani in Pakistan
- Mandenka in Senegal
- Mayan in Mexico
- Mbuti in Democratic Republic of Congo
- Miao in China
- Mongolian in China
- Mozabite in Algeria
- Naxi in China
- Northern Han in China
- Orcadian in Orkney
- Oroqen in China
- Palestinian in Israel
- Papuan in Papua New Guinea
- Pathan in Pakistan
- Pima in Mexico
- Popuan Sepik in New Guinea
- Russian in Russia
- San in Namibia
- Sardinian in Italy
- She in China
- Sindhi in Pakistan
- Surui in Brazil
- Tu in China
- Tujia in China
- Tuscan in Italy
- Uygur in China
- Xibo in China
- Yakut in Siberia
- Yi in China
- Yoruba in Nigeria


# Source
Datasets are from Stanford U, contains ~ 660,918 tag SNPs (Illumina HuHap 650k), in autosomes, chromosome X and Y, the pseudoautosomal region and mitochondrial DNA, typed across 1043 individuals from all panel populations (Li JZ et al. Science 319: 1100-4, 2008). 

## Installation

```
if (!require(devtools)) {
    install.packages("devtools")
}
devtools::install_github("Erickcufe/hgdp", ref = "erick")
```

