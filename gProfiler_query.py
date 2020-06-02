#! /usr/bin/env python3

import csv
import requests

# load gmt
def load_gmt(gmt_file):
    gmt = {}
    with open(gmt_file) as fhd:
        for line in fhd:
            line = line.strip().split('\t')
            gmt[line[0]] = set(line[2:])
    return gmt


def gProfile_query(query_gene, gmt, output_file, gene_to_symbol, output_symbol_file):
    r = requests.post(url='https://biit.cs.ut.ee/gprofiler/api/gost/profile/',
                      json={
                          'organism': 'drerio',
                          'query': list(query_gene),
                          'sources': ['GO'],
                          'no_evidences': True
                      })
    response = r.json()
    with open(output_file, 'w') as output_fhd1, \
         open(output_symbol_file, 'w') as output_fhd2:
        fhd_csv1 = csv.writer(output_fhd1)
        fhd_csv2 = csv.writer(output_fhd2)
        fhd_csv1.writerow(
            ('source', 'term_name', 'term_id', 'adjusted_p_value',
             'negative_log10_of_adjusted_p_value', 'term_size', 'query_size',
             'intersection_size', 'effective_domain_size', 'intersections'))
        fhd_csv2.writerow(
            ('source', 'term_name', 'term_id', 'adjusted_p_value',
             'negative_log10_of_adjusted_p_value', 'enrichment', 'term_size',
             'query_size', 'intersection_size', 'effective_domain_size',
             'intersections'))
        for result in response['result']:
            intersections = gmt[result['native']] & query_gene
            row = (result['source'], result['name'], result['native'],
                   result['p_value'], -np.log10(result['p_value']),
                   result['term_size'], result['query_size'],
                   result['intersection_size'],
                   result['effective_domain_size'], ','.join(intersections))
            fhd_csv1.writerow(row)
            row = (
                result['source'], result['name'], result['native'],
                result['p_value'], -np.log10(result['p_value']),
                (1.0 * result['intersection_size'] / result['query_size']) /
                (1.0 * result['term_size'] / result['effective_domain_size']),
                result['term_size'], result['query_size'],
                result['intersection_size'], result['effective_domain_size'],
                ','.join(ensGene_to_symbol[gene] if gene in
                         ensGene_to_symbol else 'NA'
                         for gene in intersections))
            fhd_csv2.writerow(row)


gmt_file_list = [
    '/mnt/Storage/home/wangwen/source/bySpecies/danRer7/gprofiler_drerio.ENSG/drerio.GO.BP.ENSG.gmt',
    '/mnt/Storage/home/wangwen/source/bySpecies/danRer7/gprofiler_drerio.ENSG/drerio.GO.CC.ENSG.gmt',
    '/mnt/Storage/home/wangwen/source/bySpecies/danRer7/gprofiler_drerio.ENSG/drerio.GO.MF.ENSG.gmt'
]

gmt = {}
for gmt_file in gmt_file_list:
    gmt.update(load_gmt(gmt_file))

ensGene_to_symbol_file = '/mnt/Storage/home/wangwen/source/bySpecies/danRer7/danRer7_ensGeneToGeneSymbol.txt'
with open(ensGene_to_symbol_file) as fhd:
    ensGene_to_symbol = {
        line.strip().split()[0]: line.strip().split()[1]
        for line in fhd
    }


gProfile_query(query_gene, gmt, output_file, ensGene_to_symbol_file, output_symbol_file)
