#!/usr/bin/env python

"""
Clustering PDB structures using BSAS algorithm.

Created by: Brian Jimenez-Garcia <brian.jimenez@bsc.es>

Barcelona Supercomputing Center
Joint BSC-IRB Research Program in Computational Biology,
Life Sciences Department
"""

import sys
import os
import itertools
import Bio.PDB
import numpy as np


def get_ca_atoms(pdb_files):
    ca_atoms = {}
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    for struct_id, pdb_file in enumerate(pdb_files):
        print "Reading CA from %s" % pdb_file
        structure = pdb_parser.get_structure(pdb_file, pdb_file)
        model = structure[0]
        ca_atoms[pdb_file] = dict()
        for chain in model:
            chain_id = chain.get_id()
            if chain_id not in ca_atoms[pdb_file].keys():
                ca_atoms[pdb_file][chain_id] = []
            for residue in chain:
                try:
                    ca_atoms[pdb_file][chain_id].append(residue['CA'])
                except:
                    ca_atoms[pdb_file][chain_id] = [residue['CA']]
    return ca_atoms


def calculate_rmsd(atoms_a, atoms_b):
    coord_a = np.array([atom.coord for atom in atoms_a])
    coord_b = np.array([atom.coord for atom in atoms_b])
    rmsd = np.sqrt(np.sum((coord_a - coord_b)**2)/len(coord_a))
    return rmsd


def clusterize(pdb_files, chain_rec, chain_lig, cutoff):
    N = len(pdb_files)
    super_imposer = Bio.PDB.Superimposer()

    clusters_found = 0
    clusters = {clusters_found : [pdb_files[0]]}

    # Read all structures CA's
    ca_atoms = get_ca_atoms(pdb_files)

    for pdb in pdb_files[1:]:
        print "pdb %s" % pdb
        in_cluster = False
        for cluster_id in clusters.keys():
            # For each cluster representative
            representative = clusters[cluster_id][0]
            ca_atoms_rec_ref = list(itertools.chain.from_iterable([ca_atoms[representative][chain] for chain in chain_rec]))
            ca_atoms_rec_target = list(itertools.chain.from_iterable([ca_atoms[pdb][chain] for chain in chain_rec]))
            # Superimposing receptor
            super_imposer.set_atoms(ca_atoms_rec_ref, ca_atoms_rec_target)
            ca_atoms_lig_ref = list(itertools.chain.from_iterable([ca_atoms[representative][chain] for chain in chain_lig]))
            ca_atoms_lig_target = list(itertools.chain.from_iterable([ca_atoms[pdb][chain] for chain in chain_lig]))
            # Calculate Ligand-RMSD
            rmsd = calculate_rmsd(ca_atoms_lig_ref, ca_atoms_lig_target)
            print "[%s and %s]" % (representative, pdb)
            print 'Superimposed RMSD is %5.3f' % (super_imposer.rms)
            print 'L-RMSD between is %5.3f' % (rmsd)
            if rmsd <= cutoff:
                clusters[cluster_id].append(pdb)
                print "PDB %s goes into cluster %d" % (pdb, cluster_id)
                in_cluster = True
                break

        if not in_cluster:
            clusters_found += 1
            clusters[clusters_found] = [pdb]
            print "New cluster %d" % clusters_found
        print

    return clusters


def write_cluster_info(clusters):
    with open('cluster.out', 'w') as output:
        for id_cluster, ids in clusters.iteritems():
            output.write("%d:%d:[%s]\n" % (id_cluster, len(ids), ','.join(ids)))


def read_pdb_list(file_name):
    pdb_list = []
    with open(file_name) as input:
        pdb_list = [line.rstrip(os.linesep) for line in input.readlines()]
    return pdb_list


def usage():
    print "Usage: %s pdb_list chain_rec chain_lig cutoff" % sys.argv[0]
    print "  Example: %s to_cluster.list A,B C 10." % sys.argv[0]


if __name__ == '__main__':

    if len(sys.argv[1:]) != 4:
        usage()
        raise SystemExit("Error: wrong command line")

    pdb_files = read_pdb_list(sys.argv[1])
    chain_rec = [chain.strip().upper() for chain in sys.argv[2].split(',')]
    chain_lig = [chain.strip().upper() for chain in sys.argv[3].split(',')]
    cutoff = float(sys.argv[4])
    print "Total number of PDB files to clusterize is %d" % len(pdb_files)
    print "Receptor chains %s" % str(chain_rec)
    print "Ligand chains %s" % str(chain_lig)
    print "Cutoff is %5.3f" % cutoff

    clusters = clusterize(pdb_files, chain_rec, chain_lig, cutoff)

    write_cluster_info(clusters)
    print "Done."
