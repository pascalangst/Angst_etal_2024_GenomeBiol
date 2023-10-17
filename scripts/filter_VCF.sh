# Analysis code for filtering of VCF file obtained from snakemake workflow available from https://github.com/pascalangst/Angst_etal_2022_MBE and described in https://doi.org/10.1093/molbev/msac264

# fix names
zcat magna_2014-2018_all_220329.vcf.gz | grep ^# | \
sed 's/M-27_spr2014/M-73_spr2014/g' | \
sed 's/K-8_spr2016/SK-58_spr2016/g' | \
sed 's/N-41b_spr2015/N-42_spr2015/g' | \
sed 's/SKW-2b_spr2018/SKN-1_spr2018/g' | \
sed 's/G-42_spr2016/G-41_spr2016/g' | \
sed 's/G-42b_spr2016/G-42_spr2016/g' | \
sed 's/SKN-2_spr2017/SKW-2_spr2017/g' | \
cat - <(zcat magna_2014-2018_all_220329.vcf.gz | grep -v ^#) | gzip > magna_2014-2018_all_220329.fix.vcf.gz

# filter samples, quality, and biallelic SNPs
zcat magna_2014-2018_all_220329.fix.vcf.gz | \
vcffilter -f "QUAL > 30 & MQ > 40 & QD > 2.0 & FS < 60" | \
vcftools --vcf - --remove-indels --max-alleles 2 --max-missing 0.001 \
--remove-indv SK-1_spr2014 \
--remove-indv G-33_spr2014 \
--remove-indv SK-44_spr2017 \
--remove-indv LON-9A_spr2018 \
--remove-indv G-10b_smr2018 \
--remove-indv HA-1_smr2014 \
--remove-indv HA-1_smr2015 \
--remove-indv HA-1_spr2014 \
--remove-indv HA-1_spr2015 \
--remove-indv N-44Peter_spr2014 \
--remove-indv SK-17A_smr2016 \
--remove-indv N-6L_smr2014 \
--remove-indv N-71_spr2015 \
--remove-indv SKW-1_spr2015 \
--remove-indv SKN-1_spr2018 \
--recode --recode-INFO-all --out magna_2014-2018_all_220329.fix.SNP_QUAL_filter

# fix GATK DP estimate
java -jar ~/programs/jvarkit/dist/vcffilterjdk.jar -e 'return new VariantContextBuilder(variant).genotypes( variant.getGenotypes().stream().map(G->{ if(!G.isHet()) return G; final int ad[]=G.getAD(); if(ad.length!=2) return G; if(ad.length==2) { final GenotypeBuilder gb = new GenotypeBuilder(G.getSampleName(),Arrays.asList(G.getAllele(0),G.getAllele(1))); final Genotype genotype1 = gb.AD(new int[] { ad[0], ad[1] }).GQ(G.getGQ()).DP( ad[0] + ad[1] ).PL(G.getPL()).make(); return genotype1; } return G; }).collect(Collectors.toList()) ).make();' magna_2014-2018_all_220329.fix.SNP_QUAL_filter.recode.vcf > magna_2014-2018_all_220329.fix.SNP_QUAL_filter.tmp
java -jar ~/programs/jvarkit/dist/vcffilterjdk.jar -e 'return new VariantContextBuilder(variant).genotypes( variant.getGenotypes().stream().map(G->{ if(!G.isHom()) return G; final int ad[]=G.getAD(); if(ad.length!=2) return G; if(ad.length==2) { final GenotypeBuilder gb = new GenotypeBuilder(G.getSampleName(),Arrays.asList(G.getAllele(0),G.getAllele(1))); final Genotype genotype1 = gb.AD(new int[] { ad[0], ad[1] }).GQ(G.getGQ()).DP( ad[0] + ad[1] ).PL(G.getPL()).make(); return genotype1; } return G; }).collect(Collectors.toList()) ).make();' magna_2014-2018_all_220329.fix.SNP_QUAL_filter.tmp > magna_2014-2018_all_220329.fix.SNP_QUAL_filter_DPadj.vcf

# mask "calls" supported by just one allele 
java -jar ~/programs/jvarkit/dist/vcffilterjdk.jar -e 'return new VariantContextBuilder(variant).genotypes( variant.getGenotypes().stream().map(G->{ if(!G.isHet()) return G; final int ad[]=G.getAD(); if(ad.length!=2) return G; final double treshold = 1; final double R= ad[0]; final double A= ad[1]; if((( A <= treshold && A > 0)) || (( R <= treshold && R > 0))) { return GenotypeBuilder.createMissing(G.getSampleName(),G.getPloidy()); } return G; }).collect(Collectors.toList()) ).make();' magna_2014-2018_all_220329.fix.SNP_QUAL_filter_DPadj.vcf > magna_2014-2018_all_220329.fix.SNP_QUAL_AD_filter.vcf
java -jar ~/programs/jvarkit/dist/vcffilterjdk.jar -e 'return new VariantContextBuilder(variant).genotypes( variant.getGenotypes().stream().map(G->{ if(!G.isHom()) return G; final int ad[]=G.getAD(); if(ad.length!=2) return G; final double treshold = 1; final double R= ad[0]; final double A= ad[1]; if((( A <= treshold && A > 0)) || (( R <= treshold && R > 0))) { return GenotypeBuilder.createMissing(G.getSampleName(),G.getPloidy()); } return G; }).collect(Collectors.toList()) ).make();' magna_2014-2018_all_220329.fix.SNP_QUAL_AD_filter.vcf > magna_2014-2018_all_220329.fix.SNP_QUAL_AD_filter_2.vcf

# max DP filter (DP < 2*mode)
mkdir samplefiltered_2014-2018_2mode

declare -A filter
filter=$(paste <(echo "filter=(") <(paste <(find ../remaining_reads/. -name *covmode) <(find ../remaining_reads/. -name *covmode -exec head {} \;) | \
sed -e 's/..\/remaining_reads\/.\/s....._201.\//["/' | sed -e 's/.bam.covmode\t/"]="/' | sed 's/$/"/' | tr "\n" " ") <(echo ")") | sed -e 's/\t//g' | \
sed 's/M-27_spr2014/M-73_spr2014/g' | \
sed 's/K-8_spr2016/SK-58_spr2016/g' | \
sed 's/N-41b_spr2015/N-42_spr2015/g' | \
sed 's/SKW-2b_spr2018/SKN-1_spr2018/g' | \
sed 's/G-42_spr2016/G-41_spr2016/g' | \
sed 's/G-42b_spr2016/G-42_spr2016/g' | \
sed 's/SKN-2_spr2017/SKW-2_spr2017/g' | \
sed 's/ )/)/')

samples=$(bcftools query -l magna_2014-2018_all_220329.fix.SNP_QUAL_AD_filter_2.vcf)

task(){
   sleep 0.5
   echo "$1"
   vcftools --vcf magna_2014-2018_all_220329.fix.SNP_QUAL_AD_filter_2.vcf --indv $1 --maxDP "$(( ${filter[$1]} * 2))" --recode --recode-INFO-all --out samplefiltered_2014-2018_2mode/$1
   bgzip -c samplefiltered_2014-2018_2mode/$1.recode.vcf > samplefiltered_2014-2018_2mode/$1.recode.vcf.gz
   bcftools index samplefiltered_2014-2018_2mode/$1.recode.vcf.gz
}

N=10

for sample in $samples
    do
   ((i=i%N)); ((i++==0)) && wait
   task "$sample" &
done

wait

bcftools merge samplefiltered_2014-2018_2mode/*vcf.gz -O z -o magna_2014-2018_all_220329.fix.SNP_QUAL_AD_MaxDP-2mode_filter.merged.vcf.gz


# min DP filter (DP > 10)
vcftools --gzvcf magna_2014-2018_all_220329.fix.SNP_QUAL_AD_MaxDP-2mode_filter.merged.vcf.gz --remove-indv N-61_smr2014 --remove-indv LON-1_smr2014 --remove-indv N-89_smr2014 --remove-indv N-71_spr2014 --minDP 10 --max-missing 0.001 --recode --recode-INFO-all --out magna_2014-2018_all_220329.fix.SNP_QUAL_AD_maxDP-2mode_filter

# make missing genotypes missing all INFOs (e.g., missing AD, DP, GQ, PL; not just GT)
cat magna_2014-2018_all_220329.fix.SNP_QUAL_AD_maxDP-2mode_filter.recode.vcf | sed -e 's~\./\.[0-9:,.]*~./.~g' > magna_2014-2018_all_220329.fix.SNP_QUAL_AD_maxDP-2mode_filter.masked.vcf
