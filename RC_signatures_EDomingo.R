# Collection of signatures related to response in CRT rectal cancer (from E. Domingo publication https://doi.org/10.1016/j.ebiom.2024.105228)
# Plots
library(xlsx)

# PART I. ----
# Ghadimi (2005)
Ghad_2005 <- c('CLMN', 'MTF1', 'FLJ12949', 'CPNE3', 'FNLB', 'FKBP1B', 'SNAPC2', 
           'MLN', 'ACAA1', 'KIAA0138', 'SMC1L1', 'PPP1R10', 'WASF2', 'RGS191P1', 
           'PAK1', 'EST', 'GPT', 'VILL', 'CDC42BPA', 
           'CCNT1', 'ITPK1', 'S164',
           'MGLL', 'HOXD9', 'ELF1', 'MLL', 'AP3D1', 'PRAC', 'C11orf13', 'MUC5B', 'KIAA0284', 'DTX2', 
           'SCNN1B', 'YF13H12', 'MYO1A', 'GUCY1B3', 'PDGFC', 
           'IL12A', 'STK18', 'SLC1A3', 'DKFZp762O076', 
           'MGC26706', 'LIV-1', 'EIF5A2', 'MBTPS2')


# Watanabe (2006) - 33 probe ids
Wat_2006 <- c('PROML1', 'PCK1', 'TDGF3', 'PPID', 'GPX2', 'MC7', 'CHIT1', 
           'HRH1', 'DUSP3', 'WDR1', 'FLNA',
           'APG12L', 'ITGA1', 'RAC2', 'TYRO3', 'COL5A2', 'LGALS1', 
           'PPP1R12A', 'COL17A1', 'NID2', 'TBCC', 'POU2AF1', 'KRT20', 'COL1A2', 
           'LUM', 'PSPHL', 'THBS2', 'COL3A1')


# Kim (2007) - 95 probe ids
Kim_2007 <- c('TGOLN2', 'ZNF84', 'FBS1', 'FLJ20674', 'YIPF2', 'PLCE1', 'GLDC', 'EPS15', 
           'SULT1A2', 'TRAF4', 'MMP14', 'RAD23B', 'ZFP276', 'MCF2L', 'MSX2', 
           'FLJ14054', 'PCSK5', 'AP2B1', 'PHF7', 'ARHGEF7', 'RPL10L', 'ULK4', 
           'RAP1A', 'GPR107', 'HSPC047', 'LOC92482', 'TFR2', 'MAX', 'P2RY13', 
           'PIP3-E', 'C1orf82', 'FRMD4A', 'TBC1D13', 'KCNH2', 'PYCARD', 'C2orf17', 
           'ASL', 'BTAF1', 'LOC283768', 'NF1', 'C7orf25', 'PDCD4', 'SEMA3C', 
           'UBR2', 'LEMD3', 'RAPGEF2', 'EZH1', 'USP14', 'HSPC152', 'C21orf45', 
           'RPS4Y1', 'TYMS', 'CPS1', 'FHL1', 'H2AFX', 'FOXO3A', 'HNRPC', 'SLC24A3', 
           'ITM2A', 'ENO1', 'TBL2', 'PDXK', 'DDX1', 'NUP62', 'IGLC2', 
           'TAF6', 'HNRPC', 'LMNB2', 'RPS2', 'BICD2', 'LOC90379', 'TPM4', 'TUBB', 
           'FLNA', 'HEMK1', 'CDC20', 'SAE1', 'C22orf18', 'EIF4A1', 'HSD3B2', 
           'MLF2', 'ELAVL1', 'HNRPR', 'ILF3', 'BMP4', 'LOC339692', 
           'HTR2C', 'FGFR4', 'SOX4')


# Park (2020) - 9 genes
Park_2020 <- c('FGFR3', 'GNA11', 'H3F3A', 'IL12A', 'IL1R1', 'IL2RB', 'NKD1', 'SGK2', 'SPRY2')


# RSS - 33 genes
RSS <- c("ADH5", "ALOX5AP", "ATF7", "C7orf50", "CACTIN-AS1", "CITED2", "CYBA", "EXOC1", 
         'FTH1P3', "GLG1", "HNRNPA0", "IL12RB1", "KLK14", "MAGED2", "MAST4", "MICU2", 
         'MTG1', "NACC1", "PHF20L1", "PPCDC", "PXT1", "RNASE4", "RND1", "RSL24D1", 
         "SNORD74", "SPR", "SURF2", "TCEAL4", "TMEM176A", "TUBA8", "USP30-AS1", 
         "YBX1P4", "ZNF585B")


# Gim (2016) - 65 genes (MI - minimal regression model) & 98 genes (TO - complete regression model) 
sheets <- openxlsx::getSheetNames('/home/melirapti/Rectal/Analysis/Main/data/OtherSignatures/Enric_Domingo_paper/Gim2016_13014_2016_623_MOESM2_ESM.xlsx')
gene.l <- lapply(sheets[1], openxlsx::read.xlsx, xlsxFile='/home/melirapti/Rectal/Analysis/Main/data/OtherSignatures/Enric_Domingo_paper/Gim2016_13014_2016_623_MOESM2_ESM.xlsx')

g1 <- gene.l[[1]][["MI.genes"]][!is.na(gene.l[[1]][["MI.genes"]])]

Gim_2016 <- unique(c(g1, gene.l[[1]][["TO.genes"]]))


# Agostini (2015) - 9 genes
Ago_2016 <- read.csv('/home/melirapti/Rectal/Analysis/Main/data/OtherSignatures/Enric_Domingo_paper/Agostini2015_t0002-10.1080_15384047.2015.1046652.csv', header = TRUE)
Ago_2016 <- Ago_2016$Gene.Symbol


# Palma (2014)
library(docxtractr)
doc <- docxtractr::read_docx("/home/melirapti/Rectal/Analysis/Main/data/OtherSignatures/Enric_Domingo_paper/Palma2014_pone.0112189.s004.docx")

## Extract the tables from the document
tables <- docx_extract_all_tbls(doc, header = TRUE)

## Adjust the index if your table is not the first one
gene_table <- tables[[1]]

Palm_2014 <- gene_table$V1 
Palm_2014 <- Palm_2014[-1] 


# Millino (2017) -241 unique genes
Mil_2017 <- unique(c('NRG1', 'ITGA2', 'KIAA1683', 'DQ786194', 'DFFB', 'THC2479309', 'THC2658419', 
                     'CASKIN2', 'KLF7', 'AK094623', 'AK023572', 'THC2503819', 'KIFC2', 'EHMT2', 
                     'ADD1', 'DST', 'AK057576', 'CBR3', 'HADHA', 'LOC399744', 'ECHDC3', 'OTUD7B', 
                     'THC2660361', 'THC2645960', 'ENST00000358431', 'AVIL', 'ALS2CL', 'ECHDC3', 
                     'CTSL2', 'A_32_P150269', 'ENST00000379447', 'IGSF9', 'ARFRP1', 'THC2559651', 
                     'MAPK7', 'LOC374491', 'ITSN1', 'C11orf56', 'DST', 'ALS2CL', 'AK092942', 
                     'KIAA2002', 'NOTCH2NL', 'RNF150', 'WWP2', 'AF289615', 'C1QTNF1', 'DUSP8', 
                     'DUSP8', 'NGEF', 'FLJ14186', 'NISCH', 'ERC1', 'FOXP4', 'NPIP', 'AL134462', 
                     'NOTCH2NL', 'DGKD', 'CTAGE4', 'NCAPH2', 'PCID2', 'WNK2', 'THC2646626', 
                     'AKAP8L', 'DISP2', 'DKFZp434F142', 'ATN1', 'AK054562', 'SLC5A11', 'PER2', 
                     'SEMA4B', 'BC092421', 'TM6SF1', 'RNF144', 'GBP3', 'AK054562', 'DENND1A', 
                     'THC2506246', 'A_24_P392661', 'THC2694422', 'CYP4F12', 'AK097700', 'AK022016', 
                     'LOC731275', 'RERE', 'MCF2L', 'THC2520478', 'LOC440353', 'NPIP', 'SLC12A7', 
                     'KIAA0284', 'FLJ45445', 'THC2572360', 'AK094623', 'MAN2C1', 'LOC401357', 
                     'AFG3L1', 'AK123446', 'ABCC2', 'MICAL1', 'A_24_P666795', 'THC2722939', 
                     'BM854107', 'FUNDC2', 'FAM92A1', 'LOC401357', 'ENST00000380310', 'VPS39', 
                     'THC2554498', 'FLJ45445', 'FLJ43339', 'EPB41L4B', 'AK091744', 'A_32_P160670', 
                     'ZFYVE28', 'AK095213', 'HGF', 'GRAMD3', 'IKBKB', 'LOC440353', 'FLJ22536', 
                     'AW172589', 'SPRY4', 'CRYBB2', 'PPARD', 'FLJ11710', 'AL833005', 'RAPGEF2', 
                     'HNF4A', 'KLF13', 'C20orf112', 'LAMA5', 'SERINC5', 'EME2', 'A_24_P600036', 
                     'RHBDF1', 'SLC44A2', 'THC2585656', 'PLXNB1', 'FLJ23556', 'A_32_P42666', 
                     'THC2663978', 'SLC37A2', 'THC2649467', 'AK026267', 'AK090397', 'PDXDC2', 
                     'AK022150', 'SH3PXD2A', 'TMPRSS3', 'FLJ43692', 'CTAGE4', 'RASEF', 'LOC285908', 
                     'TMEM16J', 'LOC339047', 'RAPGEF4', 'LOC401357', 'CACNB1', 'CLIP1', 
                     'ENST00000380946', 'STC2', 'LOC728411', 'BC041926', 'TPCN1', 'CGN', 
                     'THC2570492', 'LMO7', 'BE710245', 'ZMIZ1', 'FLJ43339', 'FAM65A', 
                     'A_23_P41824', 'A_32_P214860', 'CTAGE5', 'KRT8P10', 'AV661884', 'SRGAP1', 
                     'CYP1A2', 'NUB1', 'BU726029', 'ENG', 'ATP13A1', 'EPS8L3', 'ZNF506', 
                     'RASL11B', 'AI937300', 'PTTG1IP', 'CDC42', 'OTC', 'DNAJC10', 'SNN', 
                     'ZAK', 'THC2719403', 'BFAR', 'KBTBD8', 'MGC7036', 'CDCA2', 'THC2718727', 
                     'CDK5R1', 'ESCO2', 'ZNRD1', 'CDC23', 'KIAA1666', 'HMMR', 'BC029473', 
                     'NUP62CL', 'RPS19', 'NETO2', 'KIAA1524', 'MAD2L1', 'PBK', 'SH2D1A', 
                     'DDX19A', 'SLC25A40', 'C18orf54', 'LDHC', 'CCL21', 'NOS2A', 'PRR11', 
                     'CCDC78', 'MYCBP', 'LOC221710', 'C6orf173', 'SNHG7', 'LOC128977', 
                     'SNHG7', 'HAT1', 'TCF7', 'THC2661836', 'PPAT', 'S100A8', 'CCDC117', 
                     'GLIS3', 'PPIL6', 'BC041389', 'PDDC1', 'IL27RA', 'AK026368', 'KIAA1328', 
                     'CCL20', 'GZMK', 'IL7R', 'CLDN14', 'VSIG9', 'CMKLR1', 'NUP93', 'CCR6', 
                     'TRAM1', 'GTSE1', 'ENST00000366569', 'TMEM188', 'MYO1B', 'BI836406', 
                     'BCL2L13', 'LOC730038'))



# Agostini (2015) - 9 genes
Ago_2016 <- read.csv('/home/melirapti/Rectal/Analysis/Main/data/OtherSignatures/Enric_Domingo_paper/Agostini2015_t0002-10.1080_15384047.2015.1046652.csv', header = TRUE)
Ago_2016 <- Ago_2016$Gene.Symbol

signature.list <- list(Ago_2016 = Ago_2016,
                       g1 = g1,
                       Ghad_2005 = Ghad_2005,
                       Gim_2016 = Gim_2016,
                       Kim_2007 = Kim_2007,
                       Park_2020 = Park_2020,
                       RSS = RSS,
                       Wat_2006 = Wat_2006)
saveRDS(file = "/export/scratch/pangelin/HUG/RC_project/public.signature.list.Rds", signature.list)

# PART II. ----
library(dplyr)
library(GSVA)
library(GSEABase)

spe <- readRDS(file = "/home/melirapti/Rectal/Analysis/Main/New Analysis ML/GEO_singlecellexperiment.rds")

scaled <- as.data.frame(assay(spe, "scalelogcounts"))

df <- dplyr::select(scaled, -contains(c('GSE190826', 'GSE45404-1', 'GSE45404-2')))

mtx <- as.matrix(df)
gene_signatures <- list(Ghad_2005, Kim_2007, Wat_2006)

gsva_results <- gsva(mtx, gene_signatures, method = "gsva", kcdf = "Gaussian", mx.diff = TRUE)
