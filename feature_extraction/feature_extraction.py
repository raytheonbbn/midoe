from Bio import SeqIO
from Bio.Seq import Seq
import math

class FeatureExtractor():
    INLINE = 'http://sbols.org/v2#inline'
    REVERSE_COMPLEMENT = 'http://sbols.org/v2#reverseComplement'

    def __init__(self, json_data):
        self.json_data = json_data

    @classmethod
    def get_agents(cls, json_data):
        agents = set()

        if 'detections' in json_data:
            for detection in json_data['detections']:
                agents.add(detection['agent'])

        return list(agents)

    @classmethod
    def get_features_by_role(cls, json_data):
        role_to_features = {}

        if 'features' in json_data:
            for feature in json_data['features']:
                if 'role' in feature:
                    role = feature['role']
                else:
                    role = 'region'

                if role not in role_to_features:
                    role_to_features[role] = []

                role_to_features[role].append(feature)

        return role_to_features

    def __get_detections(self, agents=[]):
        detections = []

        if 'detections' in self.json_data:
            if len(agents) == 0:
                return self.json_data['detections']
            else:
                agents = set(agents)

                for detection in self.json_data['detections']:
                    if detection['agent'] in agents:
                        detections.append(detection)

        return detections

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

    @classmethod
    def __add_region_feature(cls, region, region_feature, region_feature_dict):
        if region_feature['id'] in region_feature_dict:
            (clash_region, clash_feature) = region_feature_dict[region_feature['id']]

            clash_elements = clash_feature['sequence'].lower()

            if (region_feature['sequence'].lower() != clash_elements
                    and Seq(region_feature['sequence']).reverse_complement().lower() != clash_elements):
                if clash_feature['id'] == region_feature['id']:
                    cls.__uniquely_add_region_feature(clash_region, clash_feature, region_feature_dict)

                cls.__uniquely_add_region_feature(region, region_feature, region_feature_dict)
        else:
            region_feature_dict[region_feature['id']] = (region, region_feature)

    def extract_detection_features(self, engineered=1, agents=[], min_length=1, max_length=0,
            fasta_paths=[], flanking=0, fixed_length=0):
        extracted_data = {'features': []}

        record_dict = {}

        for fasta_path in fasta_paths:
            with open(fasta_path, "r") as fasta_handle:
                record_dict.update(SeqIO.to_dict(SeqIO.parse(fasta_handle, "fasta")))

        feature_dict = {}

        if 'features' in self.json_data:
            for feature in self.json_data['features']:
                if feature['id'] not in feature_dict:
                    feature_dict[feature['id']] = []

                feature_dict[feature['id']].append(feature)

        region_feature_dict = {}

        for detection in self.__get_detections(agents):
            if 'contigs' in detection:
                for contig in detection['contigs']:
                    if 'regions' in contig:
                        for region in contig['regions']:
                            if ((engineered == 1 and 'engineered' in region and region['engineered'] == 1)
                                    or (engineered == 0 and ('engineered' not in region or region['engineered'] == 0))
                                    or engineered == 2):
                                if 'feature' in region and region['feature'] in feature_dict:
                                    features = feature_dict[region['feature']]
                                else:
                                    features = [{}]

                                for feature in features:
                                    if contig['source'] in record_dict:
                                        record = record_dict[contig['source']]

                                        start = region['start'] - flanking
                                        end = region['end'] + flanking

                                        if fixed_length > 0:
                                            region_length = end - start + 1

                                            if region_length < fixed_length:
                                                start = start - math.floor((fixed_length - region_length)/2)
                                                end = end + math.ceil((fixed_length - region_length)/2)
                                            elif region_length > fixed_length:
                                                start = start + math.floor((region_length - fixed_length)/2)
                                                end = end - math.ceil((region_length - fixed_length)/2)

                                        if start < 1:
                                            start = 1
                                        if end > len(record):
                                            end = len(record)

                                        if fixed_length <= 0 or end - start + 1 == fixed_length:
                                            if 'orientation' in region:
                                                if region['orientation'] == self.INLINE:
                                                    sequence = str(record[start - 1:end].seq)
                                                elif region['orientation'] == self.REVERSE_COMPLEMENT:
                                                    sequence = str(record[start - 1:end].reverse_complement().seq)
                                            else:
                                                sequence = str(record[start - 1:end].seq)
                                        else:
                                            sequence = ''
                                    elif 'sequence' in feature and flanking == 0 and (fixed_length == 0 or len(feature['sequence'] == fixed_length)):
                                        sequence = feature['sequence']
                                    else:
                                        sequence = ''

                                    if (len(sequence) > 0 and len(sequence) >= min_length
                                            and (max_length <= 0 or len(sequence) <= max_length)):
                                        if 'role' in region:
                                            role = region['role']
                                        elif 'role' in feature:
                                            role = feature['role']
                                        else:
                                            role = 'region'

                                        if 'id' in feature:
                                            extracted_id = feature['id']
                                        else:
                                            extracted_id = role

                                        extracted_feature = {
                                            'id': extracted_id,
                                            'role': role,
                                            'sequence': sequence
                                        }

                                        self.__add_region_feature(region, extracted_feature, region_feature_dict)

        for feature_ID in region_feature_dict:
            (region, region_feature) = region_feature_dict[feature_ID]

            if region_feature['id'] == feature_ID:
                extracted_data['features'].append(region_feature)

        return extracted_data
