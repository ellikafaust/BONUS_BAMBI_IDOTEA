# !/bin/bash

### NOT included here is the trimming, filtering and uniquing of the fastq files, starting here with .uni files:

#First, merge all the .uni files:
~/scripts/mergeUniq.pl uni >IB14_pruned_merged.uniq

#Clustering with cd-hit-est:
~/scripts/cd-hit-est -i mergedUniqTags.fasta -o cdh_alltags.fas -aL 1 -aS 1 -g 1 -c 0.91 -M 0 -T 0

#Merging the clustering output into alltags.ul:
~/scripts/uniq2loci.pl IB14_pruned_merged.uniq cdh_alltags.fas.clstr cdh_alltags.fas minDP=20 >alltags.ul

#Calling SNPs with haplocall_denovo, using filter settings recommended by Mischa Matz (cranked up as high_mem ran out of memory with default settings):
~/scripts/haplocall_denovo.pl alltags.ul mindp=5 ind=10 hetero=0.8 aobs=20 strbias=20 found=0.5

# OUTPUT::
# 22356494        total tags
# 7051292 pass depth cutoff 20
# 6799207 pass strand bias cutoff 5
# reading start: Thu Aug 18 16:47:18 2016
# 
# allele filters start: Thu Aug 18 19:00:52 2016
# 
# locus filters start: Thu Aug 18 22:59:57 2016
# VCF writing start: Thu Aug 18 23:44:10 2016
# Haplo-VCF writing start: Thu Aug 18 23:58:15 2016
# all done: Fri Aug 19 00:07:39 2016
# Input parameters (parameter names):
# clipping 0 bases from read ends (clip)
# must see an allele at least 20 times (aobs)
# strand bias cutoff (strbias): 20
# allele bias cutoff (abias): 10
# keep loci genotyped in at least 0.5 fraction of all samples (found)
# must see an allele in at least 10 individual(s) (ind)
# keep monomorphic tags? (mono): 0
# maximum acceptable fracton of heterozygotes at a locus (hetero): 0.8
# minimum depth to call a homozygote (mindp): 5
# 
#         file alltags.ul
# 
# -----------------
# 
# Allele filters:
# 
# 6796332 raw alleles
# 6796332 with 20 or more reads
# 5470290 in 10 or more samples
# 5288794 pass strand bias cutoff 20
# 4724303 pass allele bias cutoff 10
# 
# --------------------
# 
# Locus filters:
# 
# 3422577 total
# 3422577 remain after applying allele filters
# 3405097 have less than 0.8 fraction of heterozygotes
# 3235182         with not more than 2 alleles per individual
# 37304   genotyped in 50% of samples

# making a tab-delimited table of clone (replicate) sample pairs
nano clonepairs.tab
#edit accordingly 
# IB14BON02a      IB14BON02b
# IB14BON05a      IB14BON05b
# IB14BON11-1     IB14BON11-2
# IB14BON14-1     IB14BON14-2
# IB14ESP07A      IB14ESP07B
# IB14ESP09A      IB14ESP09B
# IB14ESP21A      IB14ESP21B
# IB14ESP30A      IB14ESP30B
# IB14FAL17A      IB14FAL17B
# IB14HAN18-1     IB14HAN18-2
# IB14HAN20-1     IB14HAN20-2
# IB14HAN23-1     IB14HAN23-2
# IB14HAN25-1     IB14HAN25-2
# IB14KIV01A      IB14KIV01B
# IB14KIV04A      IB14KIV04B
# IB14KIV05A      IB14KIV05B
# IB14KIV08A      IB14KIV08B
# IB14LET02a      IB14LET02b
# IB14LET08-1     IB14LET08-2
# IB14LET23A      IB14LET23B
# IB14PUL01A      IB14PUL01B
# IB14PUL02A      IB14PUL02B
# IB14PUL03A      IB14PUL03B
# IB14PUL05A      IB14PUL05B
# IB14SAL05a      IB14SAL05b
# IB14SAL13a      IB14SAL13b
# IB14SAL14a      IB14SAL14b
# IB14SAL15a      IB14SAL15b
# IB14SAS16-1     IB14SAS16-2
# IB14SAS17-1     IB14SAS17-2
# IB14SAS19-1     IB14SAS19-2
# IB14SAS21-1     IB14SAS21-2
# IB14STO01a      IB14STO01b
# IB14STO03a      IB14STO03b
# IB14STO07a      IB14STO07b
# IB14STO12-1     IB14STO12-2
# IB14TJA01A      IB14TJA01B
# IB14TJA02A      IB14TJA02B
# IB14TJA04A      IB14TJA04B
# IB14TJA26-1     IB14TJA26-2
# IB14VAS01A      IB14VAS01B
# IB14VAS02A      IB14VAS02B
# IB14VAS03A      IB14VAS03B
# IB14VAS05A      IB14VAS05B
# IB14BOR01a      IB14BOR01b
# IB14BOR02a      IB14BOR02b
# IB14BOR03a      IB14BOR03b
# IB14BOR04a      IB14BOR04b
# IB14DJU06a      IB14DJU06b
# IB14DJU07a      IB14DJU07b
# IB14DJU18a      IB14DJU18b
# IB14DJU24a      IB14DJU24b
# IB14HOO04a      IB14HOO04b
# IB14KIE04a      IB14KIE04b
# IB14KIE07a      IB14KIE07b
# IB14KIE08a      IB14KIE08b
# IB14KIE13a      IB14KIE13b
# IB14KUG01a      IB14KUG01b
# IB14KUG09a      IB14KUG09b
# IB14KUG12a      IB14KUG12b
# IB14NEU10a      IB14NEU10b
# IB14NEU11a      IB14NEU11b
# IB14NEU20-1     IB14NEU20-2
# IB14NEU25a      IB14NEU25b
# IB14NYN01a      IB14NYN01b
# IB14NYN02a      IB14NYN02b
# IB14NYN03a      IB14NYN03b
# IB14NYN04a      IB14NYN04b
# IB14OST07a      IB14OST07b
# IB14OST12a      IB14OST12b
# IB14OST17a      IB14OST17b
# IB14OST20a      IB14OST20b
# IB14OTT17A      IB14OTT17B
# IB14OTT18A      IB14OTT18B
# IB14OTT19A      IB14OTT19B
# IB14OTT20A      IB14OTT20B
# IB14PAK07A      IB14PAK07B
# IB14PAK08-1     IB14PAK08-2
# IB14PAK10A      IB14PAK10B
# IB14PAK11A      IB14PAK11B
# IB14PAN02-1     IB14PAN02-2
# IB14PAN11A      IB14PAN11B
# IB14PAN12A      IB14PAN12B
# IB14PAN22A      IB14PAN22B
# IB14PAN26-1     IB14PAN26-2
# IB14SAR02A      IB14SAR02B
# IB14SAR14A      IB14SAR14B
# IB14SAR28A      IB14SAR28B
# IB14SKR111-1    IB14SKR111-2
# IB14SKR113-1    IB14SKR113-2
# IB14SKR123-1    IB14SKR123-2
# IB14SKR124-1    IB14SKR124-2
# IB14SKR201a     IB14SKR201b
# IB14SYL01a      IB14SYL01b
# IB14SYL14a      IB14SYL14b
# IB14SYL15a      IB14SYL15b
# IB14SYL21a      IB14SYL21b
# IB14VEJ08a      IB14VEJ08b
# IB14VEJ09a      IB14VEJ09b
# IB14SKR123-1    IB14SKR123-2
# IB14SKR124-1    IB14SKR124-2
# IB14SKR201a     IB14SKR201b
# IB14SYL01a      IB14SYL01b
# IB14SYL14a      IB14SYL14b
# IB14SYL15a      IB14SYL15b
# IB14SYL21a      IB14SYL21b
# IB14VEJ08a      IB14VEJ08b		!!!VEJ08a has very low coverage, skipping this clonepair!!!!
# IB14VEJ09a      IB14VEJ09b
# IB14VEJ22a      IB14VEJ22b
# IB14VEJ25a      IB14VEJ25b
# IB14HEL01-1     IB14HEL01-2
# IB14HEL03-1     IB14HEL03-2
# IB14HEL12-1     IB14HEL12-2
# IB14HEL18-1     IB14HEL18-2
# IB14BJO05a      IB14BJO05b
# IB14BJO06a      IB14BJO06b
# IB14BJO07a      IB14BJO07b
# IB14BJO08a      IB14BJO08b
# IB14HAL02a      IB14HAL02b
# IB14HAL07a      IB14HAL07b
# IB14HAL21-1     IB14HAL21-2
# IB14HAL24a      IB14HAL24b
# IB14HER03a      IB14HER03b
# IB14HER04a      IB14HER04b
# IB14HER08-1     IB14HER08-2
# IB14RAU04a      IB14RAU04b
# IB14RAU05a      IB14RAU05b
# IB14RAU07a      IB14RAU07b
# IB14SOD07a      IB14SOD07b
# IB14SOD23-1     IB14SOD23-2
# IB14SOD24-1     IB14SOD24-2
# IB14SOD25-1     IB14SOD25-2


# extracting "true snps" subset (reproducible across replicates)
~/scripts/replicatesMatch.pl vcf=alltags.ul_Variants_count20_ab10_sb20_clip0.vcf replicates=clonepairs.tab polyonly=1 > vqsr.denovo.vcf
# 40475 total SNPs
# 1901 pass hets and match filters
# 1642 show non-reference alleles
# 1255 have alterantive alleles in at least 2 replicate pair(s)
# 1255 have matching heterozygotes in at least 0 replicate pair(s)
# 1255 polymorphic
# 1255 written

# same for haplotypes (NOT DOING HAPLOTYPES RIGHT NOW)
#~/scripts/replicatesMatch.pl vcf=alltags.ul_Vhap_count10_ab10_sb10_clip0.vcf replicates=clonepairs_no_granulosa.tab polyonly=1 > vqsr.vhap.vcf

# non-parametric recalibration (strand bias not used, -nosb)
# the goal is to generate a composite filter and determine the setting that gives maximum "gain"	
~/scripts/recalibrateSNPs.pl vcf=alltags.ul_Variants_count20_ab10_sb20_clip0.vcf true=vqsr.denovo.vcf -nosb -notp >denovo.recal.vcf

#OUTPUT::
# AB quantiles:
# 1      	13
# 5      	17
# 10     	22
# 15     	24
# 20     	26
# 30     	30
# 40     	33
# 50     	36
# 60     	39
# 70     	42
# 80     	45
# 90     	49
# 
# SB quantiles:
# 1      	24
# 5      	31
# 10     	37
# 15     	40
# 20     	45
# 30     	55
# 40     	61
# 50     	67
# 60     	74
# 70     	80
# 80     	84
# 90     	90
# 
# DP quantiles:
# 1      	9204   	97766
# 5      	10253  	56946
# 10     	11058  	52828
# 15     	11632  	47571
# 20     	12613  	45351
# 30     	14248  	40918
# 40     	15475  	37151
# 50     	16795  	34083
# 60     	18061  	30926
# 70     	19267  	28652
# 80     	20455  	26198
# 90     	21802  	24306
# 
# TP quantiles:
# 1      	1      	36
# 5      	2      	35
# 10     	3      	34
# 15     	3      	34
# 20     	4      	33
# 30     	6      	32
# 40     	7      	30
# 50     	9      	29
# 60     	10     	28
# 70     	12     	26
# 80     	16     	25
# 90     	18     	21
# 
# JOINT quantiles:
# 1      	0.00500
# 5      	0.02000
# 10     	0.04000
# 15     	0.06000
# 20     	0.08000
# 30     	0.12000
# 40     	0.18000
# 50     	0.25000
# 60     	0.35000
# 70     	0.42000
# 80     	0.54000
# 90     	0.64000
# 
# ------------------------
# 12.14% 	at qual <1 (11.14% gain)
# 41.87% 	at qual <5 (36.87% gain)
# 47.33% 	at qual <10 (37.33% gain)
# 52.67% 	at qual <15 (37.67% gain)
# 56.45% 	at qual <20 (36.45% gain)
# 62.39% 	at qual <30 (32.39% gain)
# ------------------------


# also recalibrating haplotype-wise calls (NOT DOING THIS RIGHT NOW)
#~/scripts/recalibrateSNPs.pl vcf=alltags.ul_Vhap_count10_ab10_sb10_clip0.vcf true=vqsr.vhap.vcf -nosb -notp >denovo.vhap.recal.vcf


###------ TO HERE THIS IS SAME AS PREVIOUS PIPELINE, BELOW THERE ARE CHANGES -----####


# We did replicates-based recalibration and the best gain was observed at 15:
# also, minimum genotype quality we want is Q15, biallelic only, present in 80% of samples:	
~/scripts/vcftools --vcf denovo.recal.vcf --minQ 15 --minGQ 15 --min-alleles 2 --max-alleles 2 --max-missing 0.8 --recode --recode-INFO-all --out denovo.filt0

# OUTPUT::
# After filtering, kept 748 out of 748 Individuals
# After filtering, kept 12641 out of a possible 40475 Sites

# Discarding loci with too many heterozygotes, which are likely lumped paralogs
# (by default, fraction of heterozygotes should not exceed maxhet=0.75)
# this step can also filter for the fraction of missing genotypes (default maxmiss=0.5)
~/scripts/hetfilter.pl vcf=denovo.filt0.recode.vcf >denovo.hetfilt.vcf
# OUTPUT::
# 12641 total loci
# 0 dropped because fraction of missing genotypes exceeded 0.5
# 21 dropped because fraction of heterozygotes exceeded 0.75
# 12620 written

# Thinning SNP dataset - leaving one snp per tag
# by default, it will leave the SNP with the highest minor allele frequency; this
# is good for ADMIXTURE or Fst analysis
# if you plan to use dadi, however, use it with the option criterion=maxDP-random
~/scripts/thinner.pl vcf=denovo.hetfilt.vcf > thinDenov.vcf

# OUTPUT::
# 12620 total loci
# 0 loci skipped because they were closer than 40
# 7954 loci selected

# evaluating accuracy (the most important one is Heterozygote discovery rate, last column) 
# based on replicates
~/scripts/repMatchStats.pl vcf=thinDenov.vcf replicates=clonepairs.tab 

# OUTPUT::
# pair   	gtyped 	match  	[ 00   	01     	11 ]   	HetMatch       	HomoHetMismatch	HetNoCall      	HetsDiscoveryRate
# IB14BON02a:IB14BON02b  	7358   	6103(82.9%)    	 [87%  	11%    	2% ]   	662    	150    	57     	0.86
# IB14BON05a:IB14BON05b  	7795   	7556(96.9%)    	 [88%  	10%    	2% ]   	747    	208    	3      	0.88
# IB14BON11-1:IB14BON11-2	7695   	7394(96.1%)    	 [87%  	10%    	3% ]   	756    	169    	10     	0.89
# IB14BON14-1:IB14BON14-2	7743   	7523(97.2%)    	 [87%  	11%    	2% ]   	842    	156    	8      	0.91
# IB14ESP07A:IB14ESP07B  	7812   	7624(97.6%)    	 [86%  	12%    	3% ]   	897    	153    	4      	0.92
# IB14ESP09A:IB14ESP09B  	7775   	7392(95.1%)    	 [87%  	11%    	2% ]   	832    	326    	4      	0.83
# IB14ESP21A:IB14ESP21B  	7688   	7294(94.9%)    	 [88%  	10%    	2% ]   	740    	185    	16     	0.88
# IB14ESP30A:IB14ESP30B  	7755   	7159(92.3%)    	 [87%  	11%    	2% ]   	792    	325    	22     	0.82
# IB14FAL17A:IB14FAL17B  	7683   	7491(97.5%)    	 [86%  	12%    	2% ]   	877    	117    	11     	0.93
# IB14HAN18-1:IB14HAN18-2	7813   	7578(97.0%)    	 [88%  	10%    	2% ]   	759    	186    	4      	0.89
# IB14HAN20-1:IB14HAN20-2	7830   	7543(96.3%)    	 [88%  	10%    	2% ]   	775    	160    	4      	0.90
# IB14HAN23-1:IB14HAN23-2	7774   	7558(97.2%)    	 [88%  	10%    	2% ]   	742    	149    	8      	0.90
# IB14HAN25-1:IB14HAN25-2	7769   	7431(95.6%)    	 [88%  	10%    	2% ]   	774    	288    	8      	0.84
# IB14KIV01A:IB14KIV01B  	7486   	7293(97.4%)    	 [88%  	10%    	2% ]   	743    	171    	9      	0.89
# IB14KIV04A:IB14KIV04B  	7830   	7540(96.3%)    	 [87%  	11%    	2% ]   	803    	189    	6      	0.89
# IB14KIV05A:IB14KIV05B  	7690   	7466(97.1%)    	 [88%  	10%    	2% ]   	762    	199    	7      	0.88
# IB14KIV08A:IB14KIV08B  	7786   	7399(95.0%)    	 [88%  	11%    	2% ]   	797    	349    	13     	0.81
# IB14LET02a:IB14LET02b  	7835   	7071(90.2%)    	 [88%  	11%    	1% ]   	777    	583    	24     	0.72
# IB14LET08-1:IB14LET08-2	7749   	7430(95.9%)    	 [88%  	10%    	2% ]   	760    	252    	7      	0.85
# IB14LET23A:IB14LET23B  	7691   	6948(90.3%)    	 [87%  	11%    	2% ]   	758    	188    	15     	0.88
# IB14PUL01A:IB14PUL01B  	7704   	7097(92.1%)    	 [87%  	10%    	2% ]   	726    	122    	17     	0.91
# IB14PUL02A:IB14PUL02B  	7709   	7120(92.4%)    	 [87%  	10%    	2% ]   	743    	143    	14     	0.90
# IB14PUL03A:IB14PUL03B  	7488   	7232(96.6%)    	 [87%  	10%    	3% ]   	732    	105    	7      	0.93
# IB14PUL05A:IB14PUL05B  	7608   	7416(97.5%)    	 [88%  	10%    	2% ]   	734    	156    	7      	0.90
# IB14SAL05a:IB14SAL05b  	7813   	7382(94.5%)    	 [88%  	10%    	2% ]   	707    	394    	8      	0.78
# IB14SAL13a:IB14SAL13b  	7870   	7454(94.7%)    	 [88%  	11%    	2% ]   	785    	410    	2      	0.79
# IB14SAL14a:IB14SAL14b  	7790   	7466(95.8%)    	 [88%  	10%    	2% ]   	734    	324    	4      	0.82
# IB14SAL15a:IB14SAL15b  	7776   	7318(94.1%)    	 [87%  	11%    	2% ]   	798    	458    	5      	0.78
# IB14SAS16-1:IB14SAS16-2	7679   	7435(96.8%)    	 [87%  	11%    	3% ]   	802    	131    	3      	0.92
# IB14SAS17-1:IB14SAS17-2	7531   	7233(96.0%)    	 [87%  	11%    	2% ]   	782    	175    	17     	0.89
# IB14SAS19-1:IB14SAS19-2	6828   	6601(96.7%)    	 [86%  	11%    	3% ]   	722    	148    	44     	0.88
# IB14SAS21-1:IB14SAS21-2	7795   	7456(95.7%)    	 [87%  	10%    	2% ]   	760    	246    	4      	0.86
# IB14STO01a:IB14STO01b  	7848   	6402(81.6%)    	 [89%  	10%    	2% ]   	614    	345    	92     	0.74
# IB14STO03a:IB14STO03b  	7692   	7467(97.1%)    	 [88%  	10%    	2% ]   	745    	219    	10     	0.87
# IB14STO07a:IB14STO07b  	7844   	7473(95.3%)    	 [88%  	10%    	2% ]   	725    	366    		0.80
# IB14STO12-1:IB14STO12-2	7784   	7464(95.9%)    	 [88%  	9%     	2% ]   	703    	309    	2      	0.82
# IB14TJA01A:IB14TJA01B  	7724   	7514(97.3%)    	 [83%  	13%    	3% ]   	1004   	203    	3      	0.91
# IB14TJA02A:IB14TJA02B  	7827   	7603(97.1%)    	 [87%  	10%    	2% ]   	765    	219    		0.87
# IB14TJA04A:IB14TJA04B  	7769   	7519(96.8%)    	 [86%  	11%    	2% ]   	837    	241    	4      	0.87
# IB14TJA26-1:IB14TJA26-2	7773   	7291(93.8%)    	 [88%  	11%    	2% ]   	788    	453    	4      	0.78
# IB14VAS01A:IB14VAS01B  	7841   	7479(95.4%)    	 [88%  	10%    	2% ]   	769    	212    	1      	0.88
# IB14VAS02A:IB14VAS02B  	7855   	6812(86.7%)    	 [87%  	10%    	2% ]   	710    	198    	38     	0.86
# IB14VAS03A:IB14VAS03B  	7835   	6978(89.1%)    	 [87%  	12%    	1% ]   	855    	764    	15     	0.69
# IB14VAS05A:IB14VAS05B  	7856   	6996(89.1%)    	 [86%  	13%    	0% ]   	920    	765    	4      	0.71
# IB14BOR01a:IB14BOR01b  	7876   	7563(96.0%)    	 [87%  	11%    	2% ]   	819    	264    	2      	0.86
# IB14BOR02a:IB14BOR02b  	7806   	7412(95.0%)    	 [88%  	10%    	2% ]   	726    	148    	4      	0.91
# IB14BOR03a:IB14BOR03b  	7805   	7598(97.3%)    	 [89%  	10%    	2% ]   	734    	167    	2      	0.90
# IB14BOR04a:IB14BOR04b  	7827   	7533(96.2%)    	 [88%  	10%    	2% ]   	762    	285    	3      	0.84
# IB14DJU06a:IB14DJU06b  	6386   	5973(93.5%)    	 [88%  	9%     	2% ]   	567    	268    	48     	0.78
# IB14DJU07a:IB14DJU07b  	7285   	6876(94.4%)    	 [88%  	10%    	2% ]   	710    	333    	32     	0.80
# IB14DJU18a:IB14DJU18b  	7641   	7382(96.6%)    	 [88%  	9%     	3% ]   	696    	152    	7      	0.90
# IB14DJU24a:IB14DJU24b  	7783   	6052(77.8%)    	 [87%  	11%    	2% ]   	687    	160    	60     	0.86
# IB14HOO04a:IB14HOO04b  	7789   	6839(87.8%)    	 [85%  	12%    	3% ]   	818    	126    	46     	0.90
# IB14KIE04a:IB14KIE04b  	7729   	7392(95.6%)    	 [87%  	11%    	2% ]   	780    	321    	7      	0.83
# IB14KIE07a:IB14KIE07b  	7556   	7177(95.0%)    	 [87%  	11%    	2% ]   	764    	363    	16     	0.80
# IB14KIE08a:IB14KIE08b  	5742   	5495(95.7%)    	 [86%  	11%    	3% ]   	608    	229    	190    	0.74
# IB14KIE13a:IB14KIE13b  	6634   	6295(94.9%)    	 [86%  	12%    	2% ]   	763    	313    	92     	0.79
# IB14KUG01a:IB14KUG01b  	7845   	6364(81.1%)    	 [95%  	5%     	1% ]   	305    	1257   	16     	0.32
# IB14KUG09a:IB14KUG09b  	7309   	5987(81.9%)    	 [95%  	4%     	1% ]   	254    	1045   	52     	0.32
# IB14KUG12a:IB14KUG12b  	7821   	6410(82.0%)    	 [95%  	4%     	1% ]   	252    	1170   	17     	0.30
# IB14NEU10a:IB14NEU10b  	7812   	7452(95.4%)    	 [87%  	10%    	2% ]   	782    	358    	2      	0.81
# IB14NEU11a:IB14NEU11b  	7182   	6621(92.2%)    	 [86%  	12%    	2% ]   	777    	358    	36     	0.80
# IB14NEU20-1:IB14NEU20-2	7802   	7463(95.7%)    	 [87%  	11%    	2% ]   	785    	331    	1      	0.83
# IB14NEU25a:IB14NEU25b  	7840   	7182(91.6%)    	 [87%  	11%    	2% ]   	793    	379    	29     	0.80
# IB14NYN01a:IB14NYN01b  	7855   	7666(97.6%)    	 [88%  	10%    	2% ]   	741    	143    	3      	0.91
# IB14NYN02a:IB14NYN02b  	7843   	7664(97.7%)    	 [88%  	9%     	2% ]   	724    	163    		0.90
# IB14NYN03a:IB14NYN03b  	7616   	7312(96.0%)    	 [88%  	10%    	2% ]   	749    	295    	16     	0.83
# IB14NYN04a:IB14NYN04b  	7823   	7590(97.0%)    	 [89%  	10%    	2% ]   	735    	217    	5      	0.87
# IB14OST07a:IB14OST07b  	6393   	6130(95.9%)    	 [87%  	11%    	2% ]   	654    	181    	101    	0.82
# IB14OST12a:IB14OST12b  	7585   	7369(97.2%)    	 [88%  	10%    	2% ]   	719    	209    	6      	0.87
# IB14OST17a:IB14OST17b  	7800   	7611(97.6%)    	 [88%  	10%    	2% ]   	755    	168    	2      	0.90
# IB14OST20a:IB14OST20b  	7789   	7458(95.8%)    	 [88%  	10%    	2% ]   	725    	206    	2      	0.87
# IB14OTT17A:IB14OTT17B  	7798   	7549(96.8%)    	 [87%  	11%    	2% ]   	815    	187    		0.90
# IB14OTT18A:IB14OTT18B  	6784   	6320(93.2%)    	 [87%  	12%    	1% ]   	737    	360    	70     	0.77
# IB14OTT19A:IB14OTT19B  	6608   	6407(97.0%)    	 [88%  	10%    	2% ]   	668    	196    	119    	0.81
# IB14OTT20A:IB14OTT20B  	6502   	6229(95.8%)    	 [88%  	10%    	2% ]   	642    	266    	109    	0.77
# IB14PAK07A:IB14PAK07B  	7628   	7378(96.7%)    	 [88%  	10%    	2% ]   	733    	248    	19     	0.85
# IB14PAK08-1:IB14PAK08-2	7201   	7004(97.3%)    	 [88%  	10%    	2% ]   	681    	189    	38     	0.86
# IB14PAK10A:IB14PAK10B  	7764   	7422(95.6%)    	 [89%  	10%    	2% ]   	707    	322    	9      	0.81
# IB14PAK11A:IB14PAK11B  	7469   	7212(96.6%)    	 [87%  	11%    	2% ]   	772    	235    	23     	0.86
# IB14PAN02-1:IB14PAN02-2	4323   	4078(94.3%)    	 [88%  	11%    	2% ]   	434    	243    	93     	0.72
# IB14PAN11A:IB14PAN11B  	5815   	1926(33.1%)    	 [89%  	10%    	2% ]   	187    	82     	398    	0.44
# IB14PAN12A:IB14PAN12B  	4860   	4625(95.2%)    	 [88%  	10%    	2% ]   	481    	225    	257    	0.67
# IB14PAN22A:IB14PAN22B  	7457   	7262(97.4%)    	 [88%  	10%    	2% ]   	750    	166    	11     	0.89
# IB14PAN26-1:IB14PAN26-2	7707   	7438(96.5%)    	 [88%  	10%    	2% ]   	754    	253    	4      	0.85
# IB14SAR02A:IB14SAR02B  	7834   	6338(80.9%)    	 [95%  	4%     	0% ]   	265    	1206   	26     	0.30
# IB14SAR14A:IB14SAR14B  	7496   	5883(78.5%)    	 [95%  	4%     	1% ]   	263    	1257   	69     	0.28
# IB14SAR28A:IB14SAR28B  	7821   	6276(80.2%)    	 [96%  	4%     	1% ]   	235    	1166   	34     	0.28
# IB14SKR111-1:IB14SKR111-2      	7825   	7638(97.6%)    	 [88%  	10%    	2% ]   	740    	182    	1      	0.89
# IB14SKR113-1:IB14SKR113-2      	7821   	7613(97.3%)    	 [88%  	10%    	2% ]   	757    	203    	2      	0.88
# IB14SKR123-1:IB14SKR123-2      	7815   	7647(97.9%)    	 [88%  	10%    	3% ]   	734    	162    	1      	0.90
# IB14SKR124-1:IB14SKR124-2      	7801   	7587(97.3%)    	 [88%  	10%    	2% ]   	770    	201    	3      	0.88
# IB14SKR201a:IB14SKR201b	4270   	3863(90.5%)    	 [87%  	11%    	2% ]   	412    	383    	440    	0.50
# IB14SYL01a:IB14SYL01b  	7395   	6953(94.0%)    	 [82%  	11%    	7% ]   	755    	441    	7      	0.77
# IB14SYL14a:IB14SYL14b  	7256   	6903(95.1%)    	 [82%  	11%    	7% ]   	771    	238    	7      	0.86
# IB14SYL15a:IB14SYL15b  	7355   	6971(94.8%)    	 [82%  	11%    	7% ]   	774    	361    	6      	0.81
# IB14SYL21a:IB14SYL21b  	7441   	7166(96.3%)    	 [82%  	11%    	7% ]   	811    	269    	1      	0.86
# IB14VEJ09a:IB14VEJ09b  	7536   	7286(96.7%)    	 [87%  	11%    	2% ]   	795    	241    	19     	0.86
# IB14SKR123-1:IB14SKR123-2      	7815   	7647(97.9%)    	 [88%  	10%    	3% ]   	734    	162    	1      	0.90
# IB14SKR124-1:IB14SKR124-2      	7801   	7587(97.3%)    	 [88%  	10%    	2% ]   	770    	201    	3      	0.88
# IB14SKR201a:IB14SKR201b	4270   	3863(90.5%)    	 [87%  	11%    	2% ]   	412    	383    	440    	0.50
# IB14SYL01a:IB14SYL01b  	7395   	6953(94.0%)    	 [82%  	11%    	7% ]   	755    	441    	7      	0.77
# IB14SYL14a:IB14SYL14b  	7256   	6903(95.1%)    	 [82%  	11%    	7% ]   	771    	238    	7      	0.86
# IB14SYL15a:IB14SYL15b  	7355   	6971(94.8%)    	 [82%  	11%    	7% ]   	774    	361    	6      	0.81
# IB14SYL21a:IB14SYL21b  	7441   	7166(96.3%)    	 [82%  	11%    	7% ]   	811    	269    	1      	0.86
# IB14VEJ09a:IB14VEJ09b  	7536   	7286(96.7%)    	 [87%  	11%    	2% ]   	795    	241    	19     	0.86
# IB14VEJ22a:IB14VEJ22b  	7776   	7567(97.3%)    	 [86%  	11%    	3% ]   	819    	206    	1      	0.89
# IB14VEJ25a:IB14VEJ25b  	7789   	7488(96.1%)    	 [87%  	10%    	2% ]   	785    	269    	3      	0.85
# IB14HEL01-1:IB14HEL01-2	7761   	7585(97.7%)    	 [88%  	10%    	2% ]   	723    	134    	3      	0.91
# IB14HEL03-1:IB14HEL03-2	7195   	7017(97.5%)    	 [87%  	10%    	2% ]   	728    	122    	25     	0.91
# IB14HEL12-1:IB14HEL12-2	7522   	7380(98.1%)    	 [87%  	10%    	2% ]   	764    	120    	14     	0.92
# IB14HEL18-1:IB14HEL18-2	7834   	7572(96.7%)    	 [88%  	10%    	2% ]   	761    	207    	3      	0.88
# IB14BJO05a:IB14BJO05b  	7097   	302(4.3%)      	 [63%  	36%    	1% ]   	110    	22     	760    	0.22
# IB14BJO06a:IB14BJO06b  	850    	803(94.5%)     	 [57%  	43%    	1% ]   	342    	46     	526    	0.54
# IB14BJO07a:IB14BJO07b  	7851   	7494(95.5%)    	 [88%  	10%    	2% ]   	754    	355    	1      	0.81
# IB14BJO08a:IB14BJO08b  	7856   	6755(86.0%)    	 [87%  	11%    	2% ]   	773    	339    	26     	0.81
# IB14HAL02a:IB14HAL02b  	7858   	6117(77.8%)    	 [86%  	12%    	2% ]   	725    	186    	53     	0.86
# IB14HAL07a:IB14HAL07b  	7854   	7444(94.8%)    	 [88%  	10%    	2% ]   	769    	374    		0.80
# IB14HAL21-1:IB14HAL21-2	3503   	3322(94.8%)    	 [85%  	13%    	2% ]   	426    	181    	423    	0.59
# IB14HAL24a:IB14HAL24b  	7767   	7383(95.1%)    	 [88%  	10%    	2% ]   	745    	384    	4      	0.79
# IB14HER03a:IB14HER03b  	7904   	7474(94.6%)    	 [81%  	19%    	0% ]   	1390   	249    	20     	0.91
# IB14HER04a:IB14HER04b  	7938   	7506(94.6%)    	 [80%  	20%    	1% ]   	1467   	257    	13     	0.92
# IB14HER08-1:IB14HER08-2	7944   	7371(92.8%)    	 [80%  	20%    	1% ]   	1449   	272    	20     	0.91
# IB14RAU04a:IB14RAU04b  	38     	27(71.1%)      	 [7%   	93%    	0% ]   	25     	10     	875    	0.05
# IB14RAU05a:IB14RAU05b  	7691   	7354(95.6%)    	 [88%  	10%    	2% ]   	714    	155    	12     	0.90
# IB14RAU07a:IB14RAU07b  	7835   	7421(94.7%)    	 [88%  	10%    	2% ]   	771    	414    	1      	0.79
# IB14SOD07a:IB14SOD07b  	7848   	6663(84.9%)    	 [88%  	11%    	2% ]   	704    	149    	37     	0.88
# IB14SOD23-1:IB14SOD23-2	7448   	6064(81.4%)    	 [87%  	11%    	2% ]   	668    	145    	54     	0.87
# IB14SOD24-1:IB14SOD24-2	6497   	5498(84.6%)    	 [86%  	11%    	2% ]   	630    	108    	67     	0.88
# IB14SOD25-1:IB14SOD25-2	5536   	4798(86.7%)    	 [86%  	13%    	1% ]   	607    	164    	161    	0.79
# 
# ------------------------
# hets called homos depth:
# lower 25%      	23
# median 		39
# upper 75%      	59

# creating a file of replicates and other poorly sequenced samples (low sites number/high-homozygosity,
# use vcftools --het to look those up) to remove.
~/scripts/vcftools --vcf thinDenov.vcf --het --out denovo.out

#OPEN denovo.out.het AND NOTE ALL INDIVIDUALS WITH LOW N_SITES, ADD THESE TO clones2remove

# creating final, thinned dataset without clones and poorly covered individuals
~/scripts/vcftools --vcf thinDenov.vcf --remove clones2remove --recode --recode-INFO-all --out denovo.filt

#OUTPUT:: After filtering, kept 622 out of 748 Individuals

# Outputting genotype matrix:
~/scripts/vcftools --vcf denovo.filt.recode.vcf --out IB14_allPOPs_SNPs --012

#And also Hardy-Weinberg statistics:
~/scripts/vcftools --vcf denovo.filt.recode.vcf --out IB14_allPOPs_SNPs --hardy

#-------------------
# Now the same recalibration and filtering for haplotype-wise calls 
# (you can skip this unless you plan coalescent analysis, or if you want to map haplotypes to a genome assembly)
~/scripts/replicatesMatch.pl vcf=alltags.ul_Vhap_count20_ab10_sb20_clip0.vcf replicates=clonepairs.tab polyonly=1 > vqsr.haplo.vcf
~/scripts/recalibrateSNPs.pl vcf=alltags.ul_Vhap_count20_ab10_sb20_clip0.vcf true=vqsr.haplo.vcf -nosb -notp >haplo.recal.vcf

# OUTPUT::
# ------------------------
# 11.25% 	at qual <1 (10.25% gain)
# 37.87% 	at qual <5 (32.87% gain)
# 45.91% 	at qual <10 (35.91% gain)
# 53.00% 	at qual <15 (38.00% gain)
# 58.90% 	at qual <20 (38.90% gain)
# 63.73% 	at qual <30 (33.73% gain)
# ------------------------

~/scripts/vcftools --vcf haplo.recal.vcf --minQ 20 --minGQ 15 --max-missing 0.8 --recode --recode-INFO-all --out denovo.vhap.filt0

# assessing quality
~/scripts/repMatchStats.pl vcf=denovo.vhap.filt0.recode.vcf replicates=clonepairs.tab 

# making final set with no replicates
~/scripts/vcftools --vcf denovo.vhap.filt0.recode.vcf --remove clones2remove --min-alleles 2 --minQ 20 --minGQ 15 --max-missing 0.8 --recode --recode-INFO-all --out denovo.vhap.filt
# After filtering, kept 9251 out of a possible 45373 Sites

# if you have a genome draft to check the percentage of truly unique loci:
# (but then again, if you have a genome draft better use reference-based GATK pipeline!)
~/scripts/vhap2fasta.pl denovo.vhap.filt.recode.vcf > dnloci.fasta
bowtie2 -L 16 -x /nobackup/data7/pierre/Ibaltica_genome_v1_0.fasta -U dnloci.fasta -f >dnloci.sam

# OUTPUT::



# For coalescent-based analysis (MIGRATE-N, IMa) - creating haplotypes dataset while 
keeping monomorphic tags ('mono=keep') 
echo 'haplocall_denovo.pl alltags.ul mindp=10 mono=keep' >hcdnm
launcher_creator.py -j hcdnm -n hcdnm -l hcdnmj
cat hcdnmj | perl -pe 's/12way 12/1way 12/' >hcdnmjob
qsub hcdnmjob

#=================================================

# see 2brad_denovo_analysis.txt for some analysis options

    Status API Training Shop Blog About 

    © 2016 GitHub, Inc. Terms Privacy Security Contact Help 

