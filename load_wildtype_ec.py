import os
import pandas as pd
import pyrosetta
from pyrosetta.toolbox.rcsb import load_fasta_from_rcsb, pose_from_rcsb, load_from_rcsb
from joblib import Parallel, delayed
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.core import pose as Pose
from difflib import SequenceMatcher
import warnings
import ast
def check_act_binding_in_canon_substring(canon_seq_common_start, canon_seq_common_end, act_bind_site):
    #Scheck that at least 3 residues lie in start to end for act_bind_site
    count = 0
    for site in act_bind_site:
        if site >= canon_seq_common_start and site <= canon_seq_common_end:
            count += 1
    if count >= 7:
        return True
    else: return False
    
def substring_indices(pose_seq, canon_seq, act_bind_site):
    #from canon_seq, find the inddices where it diverges from pose seq
    matcher = SequenceMatcher(None, pose_seq, canon_seq, autojunk=False)
    match = matcher.find_longest_match( 0, len(pose_seq), 0, len(canon_seq))
    #a match smaller than 50 automatically fails
    if match.size < 100:
        return 0, 0, 0, 0
    #INDEXED FROM 0!!!!
    pose_seq_common_start = match.a
    pose_seq_common_end = match.a + match.size - 1
    #add 1 to these to index them from 1, and properly compare active and binding site labels
    canon_seq_common_start = match.b +1
    canon_seq_common_end = match.b + match.size 
    #check that at least 3 act_bind_site lie between canon_seq_common_start and canon_seq_common_end
    if check_act_binding_in_canon_substring(canon_seq_common_start, canon_seq_common_end, act_bind_site):
        return pose_seq_common_start, pose_seq_common_end, canon_seq_common_start, canon_seq_common_end
    else: return 0, 0, 0, 0

def load_wildtype_if_possible(row, idx):
    entry = row['Entry']
    sequence = row['Sequence']
    pdbs = str(row['PDB']).split(";")
    #convert to list and sort
    active_binding_sites = ast.literal_eval(row['Active site'])
    active_binding_sites.extend(ast.literal_eval(row['Binding site']))
    active_binding_sites.sort()
    #remove any non 4-letter codes
    pdbs = [pdb for pdb in pdbs if len(pdb) == 4]
    #check if the rcsb id ever matches sequence
    for pdb in pdbs:
        #first clean pdb (without putting it in official wildtype directory)
        try:
            pose = pose_from_rcsb(pdb, ATOM= True)
            #only get chain 1
            if pose.chain_end(1) < pose.total_residue():
                pose.delete_residue_range_slow(pose.chain_end(1) + 1, pose.total_residue())
            
            pose_seq = pose.sequence()
            len_seq = len(pose_seq)
            pose_seq_common_start, pose_seq_common_end, canon_seq_common_start, canon_seq_common_end = substring_indices(pose_seq, sequence, active_binding_sites)
            #if zeroes returned, pass
            if canon_seq_common_end == 0:
                os.remove("/Users/robsonlab/Teetly/get_data_pyrosetta/" + pdb + ".clean.pdb")
            else:
                #clamp pose for common substring
                if pose_seq_common_end + 2 < len_seq:
                    pose.delete_residue_range_slow(pose_seq_common_end + 2, pose.total_residue())
                if pose_seq_common_start > 0:
                    pose.delete_residue_range_slow(1, pose_seq_common_start)
                #now save pose with the canon_seq_common_start and canon_seq_common_end in pdb name
                Pose.dump_comment_pdb("/Users/robsonlab/Teetly/tmp/" + entry + 
                                    "_" + pdb +
                                    "_" + str(canon_seq_common_start) +
                                    "_" + str(canon_seq_common_end) + ".pdb",
                                    pose)
                
                os.remove("/Users/robsonlab/Teetly/get_data_pyrosetta/" + pdb + ".clean.pdb")
                break
        except: #too lazy to check for other nuances...
            pass

if __name__ == '__main__':
    pyrosetta.init()
    csv = pd.read_csv("/Users/robsonlab/Teetly/get_data_pyrosetta/ecs_rest.csv")
    Parallel(n_jobs = -1, backend = "threading")(delayed(load_wildtype_if_possible)(row, idx) for idx, row in csv.iterrows())