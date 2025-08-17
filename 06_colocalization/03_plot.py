import pandas as pd
import numpy as np
import dask.dataframe as dd
import os
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='locuscompare plot')
    parser.add_argument('--query','-q',type=str, default='ENSSSCG00000028506')
    parser.add_argument('--snp','-s',type=str, default='6_158402719')
    parser.add_argument('--trt',type=str, default='S_FEEDCON_30T115')
    parser.add_argument('--tiss',type=str, default='Duodenum_BEST4_enterocytes')
    parser.add_argument('--outdir','-o',type=str, default='/disk191_2/yupf/GUT_RNA/04_SNPCALLING/case/YIPF1_FCR')
    parser.add_argument('--bfile','-b',type=str, default='/disk212/zz/microQTL/genomic_ref/GTE.1602')
    parser.add_argument('--nominal','-n',type=str, default='/disk191_2/yupf/GUT_RNA/04_SNPCALLING/decon_ct_eQTL/qtl/tensorQTL_nominal')
    parser.add_argument('--eqtl_pval','-ep',type=str, default='pval_nominal')
    args = parser.parse_args()

    ### global setting
    query=args.query
    target_snp=args.snp
    outdir=args.outdir
    trt=args.trt
    tiss=args.tiss
    bfile=args.bfile
    nominal_dir=args.nominal
    pval=args.eqtl_pval
    # query='ENSSSCG00000028506'
    # target_snp='6_158402719'
    # outdir='/disk191_2/yupf/GUT_RNA/04_SNPCALLING/case/YIPF1_FCR'
    # trt='S_FEEDCON_30T115'
    # tiss='Duodenum_BEST4_enterocytes'
    if "_" in target_snp:
        chr=target_snp.split('_')[0]
    if ":" in target_snp:
        chr=target_snp.split(':')[0]
    ### L2 setting
    gwas_dir='/disk201/chenzt/metaGWAS'
    # nominal_dir='/disk191_2/yupf/GUT_RNA/04_SNPCALLING/decon_ct_eQTL/qtl/tensorQTL_nominal'
    gene_ref='/disk191_2/yupf/GUT_RNA/03_EXPRESSION/stringtie/sample1_Duodenum/sample1_Duodenum.tsv'

    ### input of gwas and eqtl
    eQTL=pd.read_parquet(f'{nominal_dir}/{tiss}.cis_qtl_pairs.{chr}.parquet')
    eQTL=eQTL[eQTL['phenotype_id']==query]
    GWAS=dd.read_csv(F"{gwas_dir}/{trt}.txt.gz",sep='\t',compression='gzip')

    ### map gene loc
    genes=pd.read_csv(f'{gene_ref}',sep='\t')
    start=genes[genes['Gene ID']==query]['Start'].values[0]-1000000
    end=genes[genes['Gene ID']==query]['End'].values[0]+1000000
    chr=genes[genes['Gene ID']==query]['Reference'].values[0]
    GWAS=GWAS[(GWAS['chromosome']==f'chr{chr}') & (start < GWAS['position'].astype(int)) & (GWAS['position'].astype(int) < end)].compute()
    if genes[genes['Gene ID']==query]['Gene Name'].values[0] != '-':
        name=genes[genes['Gene ID']==query]['Gene Name'].values[0]
    else:
        name=None

    ### calculate -logp
    GWAS['logp']=-np.log10(GWAS['pvalue'])
    GWAS['rsid']=GWAS['variant_id'].str.split('_', n=2).str[:2].str.join('_')
    GWAS['chromosome']=GWAS['chromosome'].str.replace("chr","")
    eQTL['logp']=-np.log10(eQTL[pval])
    eQTL['rsid']=eQTL['variant_id'].str.replace(":","_")
    eQTL['chromosome']=eQTL['rsid'].str.split('_').str[0]
    eQTL['position']=eQTL['rsid'].str.split('_').str[1]
    eQTL['rsid']=eQTL['chromosome'].astype(str) + '_' + eQTL['position'].astype(str)
    GWAS=GWAS[['rsid','chromosome','position','logp']]
    eQTL=eQTL[['rsid','chromosome','position','logp']]

    ### calculate and map LD
    r2_out=f'{outdir}/r2'
    cmd=f'plink --r2 --bfile {bfile} --ld-window 999999 --ld-snp {target_snp} --ld-window-kb 1000 --ld-window-r2 0 --chr-set 18 --threads 20 --out {r2_out}'
    os.system(cmd)
    LD=pd.read_csv(f'{r2_out}.ld',delim_whitespace=True)
    LD=dict(zip(LD['SNP_B'],LD['R2']))
    GWAS=GWAS[GWAS['rsid'].isin(LD.keys())]
    eQTL=eQTL[eQTL['rsid'].isin(LD.keys())]
    GWAS['r2']=GWAS['rsid'].map(LD)
    eQTL['r2']=eQTL['rsid'].map(LD)
    merge=pd.merge(eQTL,GWAS,on='rsid')
    merge.to_csv(f'{outdir}/GWAS_eQTL_merge.txt',sep='\t',index=None)
    eQTL.to_csv(f'{outdir}/eQTL.txt',sep='\t',index=None)
    GWAS.to_csv(f'{outdir}/GWAS.txt',sep='\t',index=None)


    ### draw locuscompare plot
    gwas=f'{outdir}/GWAS.txt'
    eqtl=f'{outdir}/eQTL.txt'
    os.system(f'Rscript 03_locus_compare_plot.r --gwas {gwas} --eqtl {eqtl} --gene {query} --name {name} --snp {target_snp} --tiss {tiss} --trt {trt} --outdir {outdir}')