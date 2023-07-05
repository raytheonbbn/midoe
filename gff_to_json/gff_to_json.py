import logging
import gffutils
from sqlite3 import InterfaceError
from Bio import SeqIO
from Bio.Seq import Seq

import os
import argparse
import sys
import json

class DatabaseCreationError(Exception):

    def __init__(self):
        pass

    def __str__(self):
        return "Error creating database with gffutils, possibly due to use of non-standard 0-based coordinate system."

class GFFConverter():
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

    def __init__(self, fasta_path, fasta_uri):
        self.fasta_uri = fasta_uri

        with open(fasta_path, "r") as fasta_handle:
            self.__record_dict = SeqIO.to_dict(SeqIO.parse(fasta_handle, "fasta"))

    @classmethod
    def __generate_feature_ID(cls, base_ID, library_IDs):
        i = 1

        while '_'.join([base_ID, str(i)]) in library_IDs:
            i = i + 1

        return '_'.join([base_ID, str(i)])

    @classmethod
    def __uniquely_add_region_feature(cls, region, region_feature, region_feature_dict):
        feature_ID = cls.__generate_feature_ID(region_feature['id'], region_feature_dict)

        region['feature'] = feature_ID

        region_feature['id'] = feature_ID

        region_feature_dict[feature_ID] = (region, region_feature)

        logging.info('Derived Feature with id=%s', region_feature['id'])

    @classmethod
    def __resolve_region_feature(cls, region, region_feature, region_feature_dict):
        if cls.GFF_ALIAS in region and region['feature'] == region[cls.GFF_ALIAS]:
            if cls.GFF_ID in region:
                region['feature'] = region[cls.GFF_ID]

                region_feature['id'] = region[cls.GFF_ID]

                cls.__add_region_feature(region, region_feature, region_feature_dict)
            else:
                cls.__uniquely_add_region_feature(region, region_feature, region_feature_dict)

        else:
            cls.__uniquely_add_region_feature(region, region_feature, region_feature_dict)

    @classmethod
    def __add_region_feature(cls, region, region_feature, region_feature_dict):
        if region_feature['id'] in region_feature_dict:
            (clash_region, clash_feature) = region_feature_dict[region_feature['id']]

            clash_elements = clash_feature['sequence'].lower()

            if (region_feature['sequence'].lower() != clash_elements
                    and Seq(region_feature['sequence']).reverse_complement().lower() != clash_elements):
                if clash_feature['id'] == region_feature['id']:
                    logging.info('Clashed Features with id=%s', region_feature['id'])

                    cls.__resolve_region_feature(clash_region, clash_feature, region_feature_dict)
                else:
                    logging.info('Clashed Feature with id=%s', region_feature['id'])

                cls.__resolve_region_feature(region, region_feature, region_feature_dict)
        else:
            region_feature_dict[region_feature['id']] = (region, region_feature)

            logging.info('Derived Feature with id=%s', region_feature['id'])

    def convert_gff(self, gff_path, database_path, gff_uri, merge_strategy='create_unique', use_comma_multival_sep=False):
        json_data = {'detections': [], 'features': [],"class_dist": {}}

        gffdialect = gffutils.constants.dialect
        if use_comma_multival_sep:
            gffdialect['multival separator'] = ','
        else:
            gffdialect['multival separator'] = ';'

        detection_memo = {}
        contig_memo = {}

        region_feature_dict = {}
        region_class_dist={"Nat":0,"Eng":0}
        try:
            gffutils_db = gffutils.create_db(gff_path, database_path, force=True, merge_strategy=merge_strategy,
                    dialect=gffdialect)
        except InterfaceError:
            raise DatabaseCreationError()
        
        for feature in gffutils_db.all_features():
            if feature.source not in detection_memo:
                json_data['detections'].append({'agent': feature.source, 'gff': gff_uri, 'fasta': self.fasta_uri, 'contigs': []})

                detection_memo[feature.source] = len(json_data['detections']) - 1
                contig_memo[feature.source] = {}

            detection = json_data['detections'][detection_memo[feature.source]]

            if feature.seqid not in contig_memo[detection['agent']]:
                detection['contigs'].append({'source': feature.seqid, 'regions': []})

                contig_memo[detection['agent']][feature.seqid] = len(detection['contigs']) - 1

            contig = detection['contigs'][contig_memo[detection['agent']][feature.seqid]]

            region = {'start': feature.start, 'end': feature.stop}

            if feature.featuretype != '.':
                region['role'] = feature.featuretype

            if feature.strand == '+' or feature.strand == '.':
                region['orientation'] = 'http://sbols.org/v2#inline'
            elif feature.strand == '-':
                region['orientation'] = 'http://sbols.org/v2#reverseComplement'

            if self.ENGINEERED in feature.attributes:
                if feature.attributes[self.ENGINEERED][0] == 'TRUE':
                    region['engineered'] = 1
                    region_class_dist["Eng"] += 1
                elif feature.attributes[self.ENGINEERED][0] == 'FALSE':
                    region['engineered'] = 0
                    region_class_dist["Nat"] += 1
            
            if self.CONFIDENCE in feature.attributes:
                region[self.CONFIDENCE] = float(feature.attributes[self.CONFIDENCE][0])

            for key in feature.attributes:
                if key not in self.MAPPED_ATTRIBUTES and len(feature.attributes[key]) > 0:
                    if len(feature.attributes[key]) > 1:
                        region[key] = feature.attributes[key]
                    else:
                        region[key] = feature.attributes[key][0]

            contig['regions'].append(region)
            record_header = feature.seqid
            
            if  record_header in self.__record_dict:
                record = self.__record_dict[feature.seqid]

                if self.GFF_ALIAS in feature.attributes:
                    feature_ID = feature.attributes[self.GFF_ALIAS][0]
                elif self.GFF_ID in feature.attributes:
                    feature_ID = feature.attributes[self.GFF_ID][0]
                elif feature.featuretype != '.':
                    feature_ID = feature.featuretype
                else:
                    feature_ID = 'region'

                region['feature'] = feature_ID

                region_feature = {'id': feature_ID}

                if self.GFF_ALIAS in feature.attributes:
                    region_feature['source'] = feature.attributes[self.GFF_ALIAS][0]

                if feature.strand == '+' or feature.strand == '.' or feature.strand == '':
                    region_feature['sequence'] = str(record[feature.start - 1:feature.stop].seq)
                elif feature.strand == '-':
                    region_feature['sequence'] = str(record[feature.start - 1:feature.stop].reverse_complement().seq)
                    
                if feature.featuretype != '.':
                    region_feature['role'] = feature.featuretype
                    # in case engineered call was not put into gff
                    if self.ENGINEERED not in feature.attributes:
                        if region_feature['role'] in ["host"]:
                            region_class_dist["Nat"] += 1
                        else:
                            region_class_dist["Eng"] += 1

                self.__add_region_feature(region, region_feature, region_feature_dict)
            else:
                logging.warning('Failed to derive Feature from GFF feature at %s, %s. Contig %s was not found in FASTA.', str(feature.start), str(feature.stop), feature.seqid)

        for feature_ID in region_feature_dict:
            (region, region_feature) = region_feature_dict[feature_ID]

            if region_feature['id'] == feature_ID or  feature_ID == 'ref|{}|'.format(feature_ID) :
                json_data['features'].append(region_feature)
            else:
                print ("mismatching feature_id region_feature_id")
                print(feature_ID)
                print (region_feature['id'])
        json_data["class_dist"] = region_class_dist

        return json_data

def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gff_paths', nargs='+')
    parser.add_argument('-f', '--fasta_paths', nargs='+')
    parser.add_argument('-j', '--json_paths', nargs='*', default=[])
    parser.add_argument('-d', '--database_paths', nargs='*', default=[])
    parser.add_argument('-c', '--use_comma_multival_sep',  action='store_true')
    parser.add_argument('-l', '--conversion_log', nargs='?', default='')
    
    args = parser.parse_args(args)

    if len(args.conversion_log) > 0:
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)

        logging.basicConfig(level=logging.INFO, filename=args.conversion_log, filemode='w',
                            format='%(levelname)s : %(message)s')

    for i in range(0, min(len(args.gff_paths), len(args.fasta_paths))):
        print('Converting ' + args.gff_paths[i])

        gff_converter = GFFConverter(args.fasta_paths[i], os.path.basename(args.fasta_paths[i]))

        if i < len(args.database_paths):
            database_path = args.database_paths[i]
        else:
            (gff_root, gff_ext) = os.path.splitext(args.gff_paths[i])

            database_path = gff_root + '_gffutils.db'

        json_data = gff_converter.convert_gff(args.gff_paths[i], database_path, gff_uri=os.path.basename(args.gff_paths[i]),
            use_comma_multival_sep=args.use_comma_multival_sep)

        if i < len(args.json_paths):
            json_path = args.json_paths[i]
        else:
            (gff_root, gff_ext) = os.path.splitext(args.gff_paths[i])
            json_path = gff_root + '.json'



        with open(json_path, 'w') as json_handle:
            json.dump(json_data, json_handle, indent=2)

    print('Finished')

if __name__ == '__main__':
    main()