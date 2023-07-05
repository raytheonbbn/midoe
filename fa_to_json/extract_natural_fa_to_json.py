
from sqlite3 import InterfaceError
from Bio import SeqIO
from Bio.Seq import Seq

import os
from os.path import join as osp
from os.path import exists as ose
import argparse
import sys
import json


felix_exp_dir = osp(os.path.dirname(os.path.dirname(os.path.dirname(__file__))),"FELIX_experiments")
print(felix_exp_dir)
sys.path.append(felix_exp_dir)
from FELIX_data_preprocess import *
import ExptParamWrapper
from configparser import SafeConfigParser

exp_config = SafeConfigParser(os.environ)

from genomic_data_utils import DNA_Alphabet_Keras, DNA_Alphabet_Keras_with_star_end
from genomic_data_utils import DNA_Alphabet_Keras_with_N, DNA_Alphabet_Keras_Ambigous
script_name = os.path.splitext(os.path.basename(__file__))[0]

class FASTAConverter():
    ENGINEERED = 'engineer'
    CONFIDENCE = 'confidence'
    GFF_ID = 'ID'
    GFF_ROLE = 'type'
    GFF_ALIAS = 'alias'

    MAPPED_ATTRIBUTES = {
        ENGINEERED,
        CONFIDENCE,
        GFF_ROLE
    }

    def __init__(self, species_fa_dict,cfg_file):
        self.species_fa=species_fa_dict
        self.min_seq_len  =  200
        self.max_read_len =5000
        self.default_conf_file = cfg_file
        exp_config.read(self.default_conf_file)
        self.epw = ExptParamWrapper.ExptParamWrapper(exp_config,"fa_to_json")
        self.epw.max_read_len = self.max_read_len



    def get_regular_reads_from_read (self, seq):
        dna_alphabet_pos = DNA_Alphabet_Keras_Ambigous()
        print("init part of script", seq[:10])
        read_list,_ = get_regular_reads_from_seq([str(seq)], self.epw, dna_alphabet_pos)
        return read_list

    def convert_nat_fa(self):
        json_data_list = dict()
        # go through taxa
            # go through fa list
            # go through fa files
            # go through sequence snips
            # limit to 5000-bp seqs
            # filter too-short sequences
        for species, fa_file_list in self.species_fa.items():

            
            taxa_spec_region_count = 0
            json_data_list[species]=dict()
            for fa_file in fa_file_list:
                json_data = {'detections': [], 'features': []}
                region_feature_dict = dict()
                if not ose(fa_file):
                    continue
                fa_file_basename =os.path.basename(fa_file)
                fa_file_basename = os.path.splitext(fa_file_basename)[0]
                contig = dict(regions=[])
                for record in SeqIO.parse(fa_file, "fasta"):
                    print("record_id", record.id)
                    if '|' in record.id:
                        rec_id_str = record.id.strip()
                        z = rec_id_str.count('|')
                        id_f=rec_id_str.split('|')
                        gid=id_f[1]
                        if z >= 3:

                            accession_number = id_f[3]
                    else:
                        accession_number = record.id.strip('|')
                    print(accession_number)

                    if len(record.seq) < self.min_seq_len:
                        continue
                    if len(record.seq) > self.max_read_len:
                        regions = self.get_regular_reads_from_read(record.seq)
                    else:
                        regions= [record.seq]
                    taxa_spec_region_count += len(regions)
                    for region_i,region_sequence in enumerate(regions):
                        region = {'start': 0, 'end': len(region_sequence)}

                        region['engineered'] = 0
                        contig['regions'].append(region)
                        feature_ID = 'region'
                        region_name = accession_number+'_r'+str(region_i)
                        feature_ID = 'region_from_{}'.format(accession_number)
                        region['feature'] = feature_ID
                        region_feature = {'id': feature_ID}
                        region_feature['sequence'] = str(region_sequence).upper()
                        region_feature['role'] = "natural"
                        region_feature['source_fa'] = fa_file
                        region_feature['taxa'] = species
                        region_feature['access_number'] = accession_number
                        region_feature_dict[region_name] = (feature_ID,region_feature)

                for feature_ID in region_feature_dict:
                    (region, region_feature) = region_feature_dict[feature_ID]

                    #if region_feature['id'] == feature_ID:
                    json_data['features'].append(region_feature)
                    json_data["class_dist"]= {"Nat": taxa_spec_region_count, "Eng":0}

                json_data_list[species].update({fa_file_basename:json_data})

        return json_data_list

def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--fasta_to_spec_map',  type=str,dest="fasta_to_spec_map")
    parser.add_argument('-j', '--json_paths',  default=".")
    parser.add_argument('-c', '--exp_config_file', default='./DL_exp.cfg')

    args = parser.parse_args(args)
    species_fa_dict = dict()
    with open(args.fasta_to_spec_map,'r') as fh:
        for line in fh.readlines():
            w1,w2 =line.strip().split()
            current_fa_list =species_fa_dict.get(w1,[])
            current_fa_list.append(w2)
            species_fa_dict[w1] = current_fa_list

    fasta_converter = FASTAConverter(species_fa_dict, args.exp_config_file)
    json_data_dict = fasta_converter.convert_nat_fa()
    json_path_root = args.json_paths
    for spec,json_data in json_data_dict.items():
        print('bg seqs from species'.format(spec))
        for fa_source,json_content in json_data.items():
            json_path= osp(json_path_root,spec,"{}_{}_regions.json".format(spec,fa_source))
            if not ose(os.path.dirname(json_path)):
                os.makedirs(os.path.dirname(json_path))
            with open(json_path, 'w') as json_handle:
                json.dump(json_content, json_handle, indent=2)
    print('Finished')

if __name__ == '__main__':
    main()
