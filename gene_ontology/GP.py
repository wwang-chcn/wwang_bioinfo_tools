#! /usr/bin/env python3

import csv
import sys
from collections import defaultdict
from typing import Iterator, Optional

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd
import seaborn as sns
from gprofiler import GProfiler
from more_itertools import flatten
from scipy.spatial.distance import pdist, squareform


def BFS(node: int,
        adjacent_matrix: np.ndarray,
        discovered: list[int],
        ignored: Optional[set[int]] = None) -> None:
    """Perform BFS of the undiscovered portion of Graph g (representing in adjacent matrix) starting at node."""
    if ignored is None:
        ignored = set()
    node_number = adjacent_matrix.shape[0]
    level: list[int] = [node]
    while len(level) > 0:
        discovered.extend(level)
        next_level: list[int] = []
        for n in level:
            for n_ in range(node_number):
                if n_ not in discovered and n_ not in ignored and adjacent_matrix[
                        n][n_]:
                    next_level.append(n_)
        level, next_level = next_level, []


class QueryGenes():
    """Store query genes used in GoAna."""

    def __init__(self, query: dict) -> None:
        self.query_genes = list(set(flatten(query.values())))

    def __len__(self) -> int:
        """Get the number of all genes."""
        return self.query_genes.__len__()

    def get_bool_array(self, genes: list[str]) -> npt.ArrayLike:
        return np.array(
            [True if g in genes else False for g in self.query_genes])


class GoTerms():
    """Store GO terms used in GoAna."""

    def __init__(self, GO_output) -> None:
        self.go_terms = GO_output[[
            'description', 'name', 'native', 'source', 'term_size'
        ]].drop_duplicates().set_index('native')

    def __len__(self) -> int:
        """Get the number of all go terms."""
        return self.go_terms.shape[0]

    def get_bool_array(self, go_terms: list[str]) -> npt.ArrayLike:
        return np.array(
            [True if g in go_terms else False for g in self.go_terms.index])


class GoCluster():
    """GO terms cluster with eary query."""
    pass


class GoAna():
    """Performing GO analysis, output, and present the results."""

    def __init__(self,
                 query: dict,
                 organism: str,
                 ana_name: str = 'GoAna',
                 sources: list[str] = ['GO:CC', 'GO:MF', 'GO:BP']) -> None:
        self.query = query
        self.organism = organism
        self.ana_name = ana_name
        self.sources = sources
        self.get_query_gene_table()
        self.set_gProfiler()
        self.threshold = {
            'term_size': (15, 500),
            'p_value': 0.001,
            'topK': 3,
            'distance_method': 'jaccard',
            'distance': 0.3
        }

    def get_query_gene_table(self) -> None:
        """
        Set query table.
        
        Query gene table:
            Index: query genes,
            Columns: query."""
        self.query_name: list[str] = list(self.query.keys())
        self.query_genes = QueryGenes(self.query)
        self.query_genes_table: pd.DataFrame = pd.DataFrame(
            np.concatenate([
                self.query_genes.get_bool_array(genes)
                for _, genes in self.query.items()
            ]).reshape((len(self.query_name), -1)).T,
            index=self.query_genes.query_genes,
            columns=self.query_name,
        )

    def set_gProfiler(self) -> None:
        """Set g:profiler."""
        # self.gp = GProfiler(return_dataframe=True)  # Type  hint for GProfiler do not support pd.DataFrame ouput now.
        self.gp = GProfiler()

    def get_query_go_gene_tables(self) -> None:
        """
        Get the GO term gene table for each query.

        Filter GO terms using p-value and term size. GO terms will be differen in differen query.
        
        GO term gene table:
            Index:      Sorted GO term native,
            Columns:    Query genes."""
        self.go_terms = GoTerms(self.GO_output)
        self.query_go_gene_tables = {}
        for qname in self.query_name:
            temp = self.GO_output.query(
                f'p_value <= {self.threshold["p_value"]} & query == "{qname}" & {self.threshold["term_size"][0]} <= term_size <= {self.threshold["term_size"][1]}'
            ).sort_values('enrichment', ascending=False)
            self.query_go_gene_tables[qname] = pd.DataFrame(
                np.concatenate([
                    self.query_genes.get_bool_array(x)
                    for x in temp['intersections']
                ]).reshape((temp.shape[0], len(self.query_genes))),
                index=temp['native'],
                columns=self.query_genes.query_genes)

    def do_enrichment(self,
                      term_size: Optional[tuple] = None,
                      p_value: Optional[float] = None) -> None:
        """
        Perform gProfiler genrichment analysis with enrichment score calculated.
        
        GO output:
            columns: ['description', 'effective_domain_size', 'intersection_size', 'intersections', 'name', 'native', 'p_value', 'parents', 'precision', 'query', 'query_size', 'recall', 'significant', 'source', 'term_size', 'evidences']"""
        self.threshold['term_size'] = self.threshold[
            'term_size'] if term_size is None else term_size
        self.threshold['p_value'] = self.threshold[
            'p_value'] if p_value is None else p_value
        sys.stdout.write(
            'Performing GO enrichment using g:Profiler, please wait.\n')
        self.GO_output: pd.DataFrame = pd.DataFrame(
            self.gp.profile(
                organism=self.organism,
                sources=self.sources,
                query=self.query,
                all_results=True,
                no_evidences=False,
            ))
        self.GO_output['enrichment'] = self.GO_output[
            'intersection_size'] / self.GO_output['query_size'] / (
                self.GO_output['term_size'] /
                self.GO_output['effective_domain_size'])
        self.GO_output['-log10 P-value'] = -np.log10(self.GO_output['p_value'])
        sys.stdout.write('GO enrichment done,\n')
        self.get_query_go_gene_tables()

    def get_go_adjacent_matrix(self, go_gene_tables) -> dict:
        """
        Get adjacent matrix from GO term gene tables.
        
        GO term adjacent matrix: np.ndarray shape=(N,N), dtype=bool"""
        query_go_distance_condensed_matrix = {}
        query_go_adjacent_matrix = {}
        for qname in self.query_name:
            query_go_distance_condensed_matrix[qname] = pdist(
                go_gene_tables[qname],
                self.threshold['distance_method'])  # type: ignore
            query_go_adjacent_matrix[qname] = pd.DataFrame(
                squareform(query_go_distance_condensed_matrix[qname] <=
                           self.threshold['distance']))
        return query_go_adjacent_matrix

    def cluster_gen(
        self,
        adjacent_matrix: np.ndarray,
    ) -> Iterator[list[int]]:
        """Find cluster using BFS. Generate one by one.
        Input:
            node: npt.ArrayLike. Nodes should in order.
            adjacent_matrix: npt.DNArray.
        Output:
            Iteratior[list[str]]
            Output cluster one by one."""
        discovered, ignored = [], set()
        while True:
            for initiate_node in range(adjacent_matrix.shape[0]):
                if initiate_node not in discovered and initiate_node not in ignored:
                    BFS(initiate_node,
                        adjacent_matrix,
                        discovered=discovered,
                        ignored=ignored)
                    yield discovered
                    ignored = ignored | set(discovered)
                    discovered = []
            else:
                break  # all node included

    def query_go_terms(self, qname: str, go_terms: list[str]) -> list[list]:
        """Query go terms in GO output.
        
        Input:
            qname: str
            go_terms: list[str]
        Output:
            GO terms information: list[source, native, name, -log10 pvalue, enrichment, intersections, term_size, query_size, intersection_size, effective_domain_size]
        """
        return [[
            row['query'], row['source'], row['native'], row['name'],
            row['-log10 P-value'], row['enrichment'],
            ','.join(row['intersections']), row['term_size'],
            row['query_size'], row['intersection_size'],
            row['effective_domain_size']
        ] for _, row in self.GO_output.query(
            'native == @go_terms & query == @qname').sort_values(
                'enrichment', ascending=False).iterrows()]

    def save_within_query_cluster(self) -> None:
        """
        Save the GO cluster for each query to csv files."""
        for qname in self.query_name:
            with open(
                    f'{self.ana_name}_{qname}_top{self.threshold["topK"]}_cluster.csv',
                    'w') as fhd:
                f_csv = csv.writer(fhd)
                f_csv.writerow([
                    'query', 'source', 'native', 'name', '-log10 p-value',
                    'enrichment score', 'intersection genes', 'term size',
                    'quer size', 'intersection size', 'effective domain size'
                ])
                for index, cluster_index in enumerate(
                        self.within_query_cluster[qname]):
                    f_csv.writerow(['#', '-----', f'cluster {index+1}'] +
                                   ['-----'] * 7)
                    go_terms = self.query_go_gene_tables[qname].index[
                        cluster_index].to_list()
                    f_csv.writerows(self.query_go_terms(qname, go_terms))

    def merge_within_query(self,
                           topK: Optional[int] = None,
                           distance_method: Optional[str] = None,
                           distance: Optional[float] = None) -> None:
        """
        Merge the GO terms within each query category.
        
        Query cluster:
            dict[query_name: list[set[int]]]"""
        sys.stdout.write('Performing clustering within each query done.\n')
        # parameters processing
        self.threshold[
            'topK'] = self.threshold['topK'] if topK is None else topK
        self.threshold['distance_method'] = self.threshold[
            'distance_method'] if distance_method is None else distance_method
        self.threshold['distance'] = self.threshold[
            'distance'] if distance is None else distance
        # create query go adjacent matrix
        self.query_go_adjacent_matrix = self.get_go_adjacent_matrix(
            self.query_go_gene_tables)
        self.within_query_cluster = defaultdict(list)
        # output topK cluster for each query
        for qname in self.query_name:
            cluster_generator: Iterator[list[int]] = self.cluster_gen(
                self.query_go_adjacent_matrix[qname])
            i: int = 0
            try:
                while i < self.threshold['topK']:
                    cluster_index = next(cluster_generator)
                    self.within_query_cluster[qname].append(cluster_index)
                    i += 1
            except StopIteration:
                pass
        self.save_within_query_cluster()
        sys.stdout.write('Clustering within each query done.\n')

    def cluster_info_gen(self) -> Iterator[list[list]]:
        for index, go_terms in enumerate(self.between_query_cluster):
            for qname in self.query_name:
                for info in self.query_go_terms(qname, go_terms):
                    yield info

    def get_simplify_cluster_info(self) -> None:
        """
        Get the siplified GO terms information for each cluster.
        Keep the GO terms for each cluster and each query with highest enrichment.
        
        simplify cluster info: pd.DataFrame.
            columns: ['query', 'native', 'name' 'enrichment score', '-log10 P-value']"""
        simplify_cluster = []
        for index, go_terms in enumerate(self.between_query_cluster):
            one_simplify_cluster = []
            for qname in self.query_name:
                temp = []
                for info in self.query_go_terms(qname, go_terms):
                    temp.append([info[0], info[2], info[3], info[5], info[4]])
                if len(temp) == 0:
                    temp = pd.DataFrame([[qname, 'N/A', 'N/A', 0, 0]],
                                        columns=[
                                            'query', 'native', 'name',
                                            'enrichment score',
                                            '-log10 P-value'
                                        ])
                else:
                    temp = pd.DataFrame(temp,
                                        columns=[
                                            'query', 'native', 'name',
                                            'enrichment score',
                                            '-log10 P-value'
                                        ]).sort_values('enrichment score',
                                                       ascending=False).head(1)
                one_simplify_cluster.append(temp)
            one_simplify_cluster = pd.concat(one_simplify_cluster,
                                             ignore_index=True).sort_values(
                                                 'enrichment score',
                                                 ascending=False)
            one_simplify_cluster['native'] = [one_simplify_cluster.iloc[0, 1]
                                              ] * len(self.query_name)
            one_simplify_cluster['name'] = [one_simplify_cluster.iloc[0, 2]
                                            ] * len(self.query_name)
            simplify_cluster.append(one_simplify_cluster)
        self.simplify_cluster = pd.concat(simplify_cluster, ignore_index=True)

    def cluster_gos_info_gen(
            self, cluster_info: pd.DataFrame) -> Iterator[pd.DataFrame]:
        flag = False
        cluster = []
        for index, row in cluster_info.iterrows():
            if flag and row['query'] == self.query_name[
                    0]:  # output when encounter the first query name in first times
                yield pd.DataFrame(cluster)
                cluster = []
                flag = False
            elif row['query'] == self.query_name[
                    -1]:  # reset flag when encounter last query name
                flag = True
            cluster.append(row)
        if cluster:
            yield pd.DataFrame(cluster)

    def qname_GOs_info_gen(
            self, cluster_gos_info: pd.DataFrame
    ) -> Iterator[tuple[str, list[list]]]:
        qname_GOs_info, qname = [], None
        for index, row in cluster_gos_info.iterrows():
            if qname and row['query'] != qname:
                yield qname, qname_GOs_info
                qname_GOs_info = []
                qname = row['query']
            qname_GOs_info.append(row.to_list())
        if qname:
            yield qname, qname_GOs_info

    def save_cluster(self) -> None:
        with open(f'{self.ana_name}_cluster_info.csv', 'w') as fhd:
            f_csv = csv.writer(fhd)
            f_csv.writerow(self.cluster_info.columns)
            for index, cluster_gos_info in enumerate(
                    self.cluster_gos_info_gen(self.cluster_info)):
                f_csv.writerow(['#', '-----', '-----', f'cluster {index+1}'] +
                               ['-----'] * (self.cluster_info.shape[1] - 4))
                for qname, qname_GOs_info in self.qname_GOs_info_gen(
                        cluster_gos_info):
                    f_csv.writerow(
                        ['#', '-----', '-----', f'cluster {index+1}: {qname}'
                         ] + ['-----'] * (self.cluster_info.shape[1] - 4))
                    f_csv.writerows(qname_GOs_info)

    def within_query_cluster_index_to_go_native(
            self, cluster_index: list[int]) -> list[str]:
        """
        Convert index of within query cluster to GO terms native."""
        # Step1. within query cluster index -> GO terms index
        go_terms_index = np.unique(
            np.concatenate([
                np.where(self.query_cluster_go[c_index])
                for c_index in cluster_index
            ],
                           axis=None))
        return self.go_terms.go_terms.index[go_terms_index].to_list()

    def merge_between_query(self) -> None:
        """Merge the GO terms between query categories.

        between_query_cluster: list[pd.DataFrame[query_category_name,['name', '-log10 P-value', 'enrichment score', 'query']]]
        query cluster go table: dict[cluster * GO_term]"""
        # query cluster go table
        self.query_cluster_go = np.concatenate([
            np.concatenate([
                self.go_terms.get_bool_array(self.query_go_gene_tables[qname].
                                             index[go_terms_indexes].to_list())
                for go_terms_indexes in self.within_query_cluster[qname]
            ]) for qname in
            self.query_name  # (cluster in all querys) * len(go_terms)
        ]).reshape((-1, len(self.go_terms)))
        # query cluster adjacent matrix
        self.query_cluster_adjacent_matrix: np.ndarray = np.zeros(
            shape=(self.query_cluster_go.shape[0],
                   self.query_cluster_go.shape[0]),
            dtype=bool)
        for i in range(self.query_cluster_go.shape[0]):
            for j in range(self.query_cluster_go.shape[0]):
                if len(
                        set(self.go_terms.go_terms.index[
                            self.query_cluster_go[i]].to_list())
                        & set(self.go_terms.go_terms.index[
                            self.query_cluster_go[j]].to_list())) > 0:
                    self.query_cluster_adjacent_matrix[i][j] = True
        # between query cluster
        self.between_query_cluster: list[list[str]] = []
        discovered, ignored = [], set()
        for initiate_index in range(self.query_cluster_go.shape[0]):
            if initiate_index not in discovered and initiate_index not in ignored:
                BFS(initiate_index, self.query_cluster_adjacent_matrix,
                    discovered, ignored)
                self.between_query_cluster.append(
                    self.within_query_cluster_index_to_go_native(discovered))
                ignored = ignored | set(discovered)
                discovered = []
        # cluster information
        self.cluster_info = pd.DataFrame(
            self.cluster_info_gen(),
            columns=[
                'query', 'source', 'native', 'name', '-log10 p-value',
                'enrichment score', 'intersection genes', 'term size',
                'quer size', 'intersection size', 'effective domain size'
            ])
        self.save_cluster()
        self.get_simplify_cluster_info()

    def plot(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
        sns.scatterplot(
            x='query',
            y='name',
            size='enrichment score',
            hue='-log10 P-value',
            data=self.simplify_cluster,
            palette='YlGnBu',
            sizes=(0, 700),
            edgecolor='k',
            ax=ax,
        )
        ax.legend(bbox_to_anchor=(1, 0), loc='lower left')
        ax.set_xlabel('Gene groups')
        ax.set_ylabel('GO term name')

    def perform_all(self, ax=None) -> None:
        self.do_enrichment()
        self.merge_within_query()
        self.merge_between_query()
        self.plot(ax)


__all__ = ['GoAna']