import torch
from torch.utils.data import Dataset, DataLoader
import numpy as np
import pandas as pd
import os
from struc_feat import create_protein_graph, ego_label_set  , load_pdb
import os
from torch_geometric.data import Data
import re
import ast
from multiprocessing import Pool
from joblib import Parallel, delayed
from Bio.PDB import PDBParser


#TODO ADD EC NUMBER AS LABEL TOO
def generate_geometric_object(row):
    id = row['Entry']
    r = re.compile(f'{id}*')
    file_name = list(filter(r.match, file_list))
    try:
        if file_name and not os.path.exists(f'/Users/robsonlab/Teetly/wildtype_pts/{file_name[0]}.pt'):
            file_name = file_name[0]
            #get the graph object
            pdb_path = f'/Users/robsonlab/Teetly/wildtype_pdbs/{file_name}'
            
            #get proper active binding site range
            active_binding_sites = ast.literal_eval(row['Active site'])
            active_binding_sites.extend(ast.literal_eval(row['Binding site']))
            active_binding_sites.sort()
            start_site = int(file_name.split('_')[2])
            end_site = int(file_name.split('_')[3].split('.')[0])
            active_binding_sites = [site for site in active_binding_sites if site >= start_site and site <= end_site]

            #parse structure
            parser = PDBParser()
            protein = parser.get_structure(id, pdb_path)
            protein = protein[0]
            #sanity check
            residues = list(protein.get_residues())
            assert len(residues) == end_site - start_site + 1
            node_labels, coords, lrfs, residue_one_hot = load_pdb(protein)
            
            #get attributes so far for geometric object
            residue_one_hot = torch.tensor(residue_one_hot, dtype = torch.float)
            pos = torch.tensor(coords, dtype=torch.float)
            angle_geom = torch.tensor(lrfs, dtype=torch.float)
            
            #make graph to set edges
            graph  = create_protein_graph(active_binding_sites, protein, start_site)
            #get data.y
            y =  torch.tensor([att["y"] for node, att in graph.nodes(data=True)], dtype=torch.float)
            
            #get the data.y labels  and the one-hot encoding of each node
            label_graphs = ego_label_set(graph, active_binding_sites, start_site)
            edge_index = torch.LongTensor(list(graph.edges)).t().contiguous()
            
            #ec nums
            ecs = ast.literal_eval(row['EC_Shortened'])
            ec_vals = torch.zeros(8, dtype=torch.float)
            if len(ecs) > 0:
                for ec in ecs:
                    ec_vals[ec] = 1
            #properties
            data = Data(
                protein_id = id,
                ec_number = ec_vals,
                edge_index = edge_index,  
                label_graphs = label_graphs, 
                x = residue_one_hot, 
                pos = pos,
                angle_geom = angle_geom,
                y = y,
            )
            torch.save(data, out_dir + "/" + file_name + ".pt")
    except:
        print(f"Could not get pt for {id}")

if __name__ == '__main__':
    #convert tsv to csv
    csv = pd.read_csv("/Users/robsonlab/Teetly/get_data_pyrosetta/ec_1_4_start_truncated.csv")
    out_dir = '/Users/robsonlab/Teetly/wildtype_pts'
    file_list = list(os.listdir("/Users/robsonlab/Teetly/wildtype_pdbs"))
    
    #need to streamline these thre alues
    Parallel(n_jobs=-1, backend="threading")(delayed(generate_geometric_object)(row) for idx, row in
                                             csv.iterrows())