import sys
import os
from os import path,listdir
import subprocess
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import re 
import gzip
import networkx as nx
from networkx.algorithms.traversal.depth_first_search import dfs_tree
from networkx.drawing.nx_agraph import graphviz_layout
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

##################################
### General methods
##################################

def Call(s):
    subprocess.call(s, shell=True, executable='/bin/bash')
    
    
def bed_sort(input_file):
    call = "sort -k1,1 -k2,2n -o '%s' '%s'"%(input_file,input_file)    
    Call(call)
    
    
def gtf_sort(input_file):
    call = "sort -k1,1 -k4,4n -o '%s' '%s'"%(input_file,input_file)    
    Call(call)
        
        
def microparser(line):
    if len(line.strip())==0:
        return None
    X = line.strip().split('\t')        
    r = {a.split(' ')[0]:a.split(' ')[1].replace('"','') for a in [x.strip() for x in X[8].split(';') if len(x.strip())>0]}
    r.update({'chr':X[0],'start':X[3],'end':X[4],'strand':X[6],'feature':X[2]})
    return r


def parse_meta(text):
    assembly = re.findall('human genome \((.*)\),',text)[0].strip()
    version = re.findall('version (.*) \(',text)[0].strip()
    ensembl = re.findall('Ensembl (.*)\)',text)[0].strip()
    _date = re.findall('date: (.*)',text)[0].strip()
    return {'gencode':version,'assembly':assembly,'ensembl':ensembl,'date':_date}


##################################
### Main class definition
##################################

class GHTracker:
    def __init__(self,download_path=None,index_path=None,LATEST_VERSION=37):
        self.download_path = download_path
        self.index_path = index_path
        self.LATEST_VERSION = LATEST_VERSION
    
    ##################################
    ### Class properties
    ##################################
    @property
    def download_path(self):        
        return self._download_path
    
    @download_path.setter
    def download_path(self,download_path):
        if download_path is None:
            d = path.abspath(os.getcwd())            
        else:
            d = path.abspath(download_path)
        self._download_path = d
        self.raw_path = path.join(self.download_path,'raw')
        self.gene_path = path.join(self.download_path,'gene')
        self.meta_path = path.join(self.download_path,'meta')
#         for folder in [self.download_path, self.raw_path, self.gene_path, self.meta_path]:
#             if not path.exists(folder):
#                 os.makedirs()        
        
    @property
    def index_path(self):        
        return self._index_path
    
    @index_path.setter
    def index_path(self,index_path):
        if index_path is None:
            d = path.abspath(os.getcwd())
        else:
            d = path.abspath(index_path)
        #Only reload data structures if new path is different
        if (hasattr(self, '_index_path')) and (self._index_path == d):
            return
        self._index_path = d        
        self.ensembl_eq = path.join(self.index_path,'ensembl_equivalence.csv')
        self.gene_index = path.join(self.index_path,'gene_index.csv.gz')
        self.gene_graph = path.join(self.index_path,'gene_network.gml')
        self.df_index = None
        self.df_index_last = None
        self.df_ensembl = None
        self.G = None
        self.GR = None
        
    def load_index(self,index_path=None):
        if index_path is not None:
            self.index_path = index_path
        #TODO: CHECK IF INDEX FILES EXIST
        # Load if is not loaded before or index_path changed
        if self.df_index is None:
            df = pd.read_csv(self.gene_index)
            df = df[['gene_shortid','gene_id','gene_name','hgnc_id','gene_type','coord','gencode']]
            df2 = pd.read_csv(self.ensembl_eq)
            df = pd.merge(df,df2,left_on='gencode',right_on='gencode')
            dfx = df.sort_values(by=['gene_shortid','gencode'],ascending=[True,False]).drop_duplicates(subset=['gene_shortid'],keep='first')
            self.df_index = df
            self.df_index_last = dfx
            self.df_ensembl = df2
        print('Entries:',len(self.df_index))
        
    def load_graph(self,index_path=None):
        if index_path is not None:
            self.index_path = index_path
        #TODO: CHECK IF INDEX FILES EXIST
        # Load if is not loaded before or index_path changed
        if self.G is None:
            G = nx.read_gml(self.gene_graph)
            GR = G.reverse()
            self.G = G
            self.GR = GR
        print('Nodes:',len(self.G.nodes()),'Edges:',len(self.G.edges()))
        
    ##################################
    ### Query methods
    ##################################
    def query_by_geneid_list(self,gene_list,output_tsv=None):
        DF = self.df_index_last
        LATEST_VERSION = self.LATEST_VERSION
        gene_set = set(gene_list)
        dfx = DF[(DF['gene_id'].isin(gene_set))|(DF['gene_shortid'].isin(gene_set))].copy()    
        dfx['deprecated'] = dfx['gencode'].apply(lambda x:int(x!=LATEST_VERSION))
        dfx = dfx.reset_index(drop=True)
        if output_tsv is not None:
            dfx.to_csv(output_tsv,sep='\t',index=None)
        return dfx

    def query_geneid_history(self,gene_id,output_tsv=None):
        DF = self.df_index
        dfx = DF[(DF['gene_id']==gene_id)|(DF['gene_shortid']==gene_id)].copy()    
        dfx = dfx.sort_values(by=['gene_shortid','gencode'],ascending=[True,False])    
        dfx = dfx.reset_index(drop=True)
        if output_tsv is not None:
            dfx.to_csv(output_tsv,sep='\t',index=None)
        return dfx

    def query_genename_history(self,gene_name,output_tsv=None):
        DF = self.df_index
        dfx = DF[(DF['gene_name']==gene_name)].copy()
        short_ids = set(dfx['gene_shortid'].unique())
        dfx = DF[DF['gene_shortid'].isin(short_ids)]
        dfx = dfx.sort_values(by=['gene_shortid','gencode'],ascending=[True,False])
        dfx = dfx.reset_index(drop=True)
        if output_tsv is not None:
            dfx.to_csv(output_tsv,sep='\t',index=None)
        return dfx

    def putative_mappings(self,gene_list):
        G = self.G
        R = {}
        for gene in gene_list:
            g1 = dfs_tree(G,gene)        
            latest = None
            for x in g1.nodes():
                a = G.nodes()[x]
                #print(x,a)
                if a['status']=='latest':        
                    if latest is None:
                        latest=[x]
                    else:
                        latest.append(x)        
            #R[gene] = '' if latest is None else ','.join(latest)
            R[gene] = set() if latest is None else set(latest)
        #dfx = pd.DataFrame.from_dict(R,orient='index').reset_index()
        #dfx.columns = ['query','mappings']
        #return dfx
        return R
    
    ##################################
    ### Download methods
    ##################################
    
    
    ##################################
    ### Building methods
    ##################################