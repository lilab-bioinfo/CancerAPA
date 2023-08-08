
import sys
import os
import gzip
from collections import defaultdict


def split_qtl(tissue_name):
    qtl_file = "~/Data/GTEx_v8/aQTL/hg38/"+tissue_name+".cis_aqtl_all_approachb.txt.gz"
    line_count=0
    gene_dict={}
    gene_dict_gcta=defaultdict(list)
    gene_dict_caviar=defaultdict(list)
    snp_dict={}
    f = gzip.open(qtl_file,'rt')
    for qtl_line in f:
        if line_count>0:
            qtl_lines = qtl_line.rstrip().split("\t")
            if float(qtl_lines[5])<0.05:
                    gene_name = qtl_lines[1]
                    snp_infos = qtl_lines[0]
                    chr = snp_infos.split("_")[0]
                    pos = snp_infos.split("_")[1]
                    allel1 = snp_infos.split("_")[2]
                    allel2 = snp_infos.split("_")[3]
                    if gene_name not in gene_dict:
                            gcta_id=chr+"_"+pos+"_"+allel1+"_"+allel2+"_b38"
                            gene_dict_gcta[gene_name]=[gcta_id]
                            gene_dict_caviar[gene_name]=[[gcta_id,qtl_lines[3]]]
                            gene_dict[gene_name]=1
                    elif gene_name in gene_dict:
                            gcta_id=chr+"_"+pos+"_"+allel1+"_"+allel2+"_b38"
                            gene_dict_gcta[gene_name].append(gcta_id)
                            gene_dict_caviar[gene_name].append([gcta_id,qtl_lines[3]])
                            gene_dict[gene_name]=gene_dict[gene_name]+1
        else:
            line_count=1

    d = dict((k, tuple(v)) for k, v in gene_dict_gcta.iteritems())
    caviar = dict((k, tuple(v)) for k, v in gene_dict_caviar.iteritems())
    print(caviar)
    for current_gene in d:
        gene_directory="${WORK_DIR}/output/"+tissue_name+"/"+current_gene
        if not os.path.exists(gene_directory):
            os.makedirs(gene_directory)
        output_file_caviar=gene_directory+"/"+tissue_name+"_3aQTL_all_ttest.txt"
        output_file_caviar_file=open(output_file_caviar, "w")

        output_file_gcta=gene_directory+"/"+tissue_name+"_3aQTL_all.snp"
        output_file_gcta_file=open(output_file_gcta, "w")
        for qtl_id in caviar[current_gene]:
            output_file_caviar_file.writelines(qtl_id[0] + "\t" + qtl_id[1] + "\n")

        for qtl_id in d[current_gene]:
            output_file_gcta_file.writelines(qtl_id+"\n")


if __name__ == '__main__':
    tissue_name = sys.argv[1]
    split_qtl(tissue_name)
