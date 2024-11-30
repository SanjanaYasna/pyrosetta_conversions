import os
import pandas as pd
import pyrosetta
from pyrosetta.toolbox.rcsb import load_fasta_from_rcsb, pose_from_rcsb, load_from_rcsb
from joblib import Parallel, delayed
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.core import pose as Pose
from difflib import SequenceMatcher
import warnings
warnings.filterwarnings("ignore")

def check_act_binding_in_canon_substring(pose_seq, canon_seq, act_bind_site):
    pass

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
    return pose_seq_common_start, pose_seq_common_end, canon_seq_common_start, canon_seq_common_end

def load_wildtype_if_possible(row, idx):
    entry = row['Entry']
    sequence = row['Sequence']
    pdbs = str(row['PDB']).split(";")
    #remove any non 4-letter codes
    pdbs = [pdb for pdb in pdbs if len(pdb) == 4]
    #check if the rcsb id ever matches sequence
    for pdb in pdbs:
        #first clean pdb (without putting it in official wildtype directory)
        try:
            pose = pose_from_rcsb(pdb, ATOM= True)
        except: 
            print("Error in pose_from_rcsb for pdb: ", pdb)
            continue
        #only get chain 1
        if pose.chain_end(1) < pose.total_residue():
            pose.delete_residue_range_slow(pose.chain_end(1) + 1, pose.total_residue())
        pose_seq = pose.sequence()
        len_seq = len(pose_seq)
        pose_seq_common_start, pose_seq_common_end, canon_seq_common_start, canon_seq_common_end = substring_indices(pose_seq, sequence, None)
        #if zeroes returned, pass
        if canon_seq_common_end == 0:
            os.remove("/Users/robsonlab/Teetly/get_data_pyrosetta/" + pdb + ".clean.pdb")
            continue
        else:
            try:
                #clamp pose for common substring
                if pose_seq_common_end + 2 < len_seq:
                    pose.delete_residue_range_slow(pose_seq_common_end + 2, pose.total_residue())
                if pose_seq_common_start > 0:
                    pose.delete_residue_range_slow(1, pose_seq_common_start)
                #now save pose with the canon_seq_common_start and canon_seq_common_end in pdb name
                Pose.dump_comment_pdb("/Users/robsonlab/Teetly/wildtype_pdbs/" + entry + 
                                    "_" + pdb +
                                    "_" + str(canon_seq_common_start) +
                                    "_" + str(canon_seq_common_end) + ".pdb",
                                    pose)
                os.remove("/Users/robsonlab/Teetly/get_data_pyrosetta/" + pdb + ".clean.pdb")
                break
            except: #too lazy to check for other nuances...
                pass
        
        # #check if the sequence is equal to canon sequence
        # if sequence == pose.sequence():
        #     pose_from_rcsb(pdb, entry = entry, out_file = "/Users/robsonlab/Teetly/wildtype_pdbs/", ATOM=True)
        #     #remove the clean.pdb in current directory
        #     os.remove("/Users/robsonlab/Teetly/get_data_pyrosetta/" + pdb + ".clean.pdb")
        #     print("\n\nSEQUENCE MATCHED\n\n")
        #     #break from loop
        #     break
        # else:
        #     os.remove("/Users/robsonlab/Teetly/get_data_pyrosetta/" + pdb + ".clean.pdb")
        #    # os.remove("/Users/robsonlab/Teetly/get_data_pyrosetta/" + pdb + ".pdb")

if __name__ == '__main__':
    pyrosetta.init()
    csv = pd.read_csv("/Users/robsonlab/Teetly/get_data_pyrosetta/non_enz")
    Parallel(n_jobs = -1, backend = "threading")(delayed(load_wildtype_if_possible)(row, idx) for idx, row in csv.iterrows())