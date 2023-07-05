import argparse
import sys
import os

import json
from json_source_map import calculate
from jsonschema import validate

# Might add this exception back in once there is support for a log file can report warnings and not just errors
# There will almost certainly be cases where the coordinates length does not match the feature sequence length
# Will probably want a threshold for the difference above which a warning is given

# class CoordinatesLengthMismatch(Exception):

#     def __init__(self, start_key, end_key, start_line, end_line, seq_line):
#         self.start_key = start_key
#         self.end_key = end_key
#         self.start_line = start_line
#         self.end_line = end_line
#         self.seq_line = seq_line

#     def __str__(self):
#         return "The length specified by the '{sk}' property on line {sl} and the '{ek}' property on line {el} does not equal the length of the 'sequence' property on line {ql}.".format(sk=self.start_key,
#                                                                                                                                                                                          sl=str(self.start_line),
#                                                                                                                                                                                          ek=self.end_key,
#                                                                                                                                                                                          el=str(self.end_line),
#            

class InconsistentHostTaxa(Exception):

    def __init__(self, host_line, sub_host_line):
        self.host_line = host_line
        self.sub_host_line = sub_host_line

    def __str__(self):
        return "The 'hostTaxa' property on line {sl} contains values that do not appear in the 'hostTaxa' property on line {hl}.".format(sl=self.sub_host_line,
                                                                                                                                         hl=self.host_line)

class CoordinatesMisordered(Exception):

    def __init__(self, start_key, end_key, start_line, end_line):
        self.start_key = start_key
        self.end_key = end_key
        self.start_line = start_line
        self.end_line = end_line

    def __str__(self):
        return "The value of the '{sk}' property on line {sl} is greater than the value of the '{ek}' property on line {el}.".format(sk=self.start_key,
                                                                                                                                     sl=str(self.start_line),
                                                                                                                                     ek=self.end_key,
                                                                                                                                     el=str(self.end_line))

class CoordinatesMissingSubject(Exception):

    def __init__(self, object_class, subject_key, region_line, coord_keys):
        self.object_class = object_class
        self.subject_key = subject_key
        self.region_line = region_line
        self.coord_keys = coord_keys

    def __str__(self):
        if len(self.coord_keys) > 1:
            coord_key_str = ', '.join(["'" + coord_key + "'" for coord_key in self.coord_keys[:-1]])
        else:
            coord_key_str = "'" + self.coord_keys[0] + "'"

        if len(self.coord_keys) > 2:
            coord_key_str = ', and '.join([coord_key_str, "'" + self.coord_keys[-1] + "'"])
        elif len(self.coord_keys) > 1:
            coord_key_str = ' and '.join([coord_key_str, "'" + self.coord_keys[-1] + "'"])

        return "The {oc} object on line {rl} has properties {ck} but is missing a '{sk}' property.".format(oc=self.object_class,
                                                                                                           rl=str(self.region_line),
                                                                                                           ck=coord_key_str,
                                                                                                           sk=self.subject_key)

class FeatureReferenceError(Exception):

    def __init__(self, feature_ref, property_key, feature_line):
        self.feature_ref = feature_ref
        self.property_key = property_key
        self.feature_line = feature_line

    def __str__(self):
        return "The value '{fr}' of the '{pk}' property on line {fl} does not identify a feature object in the root list of feature objects.".format(fr=self.feature_ref,
                                                                                                                                                     pk=self.property_key,
                                                                                                                                                     fl=str(self.feature_line))

class FeatureIDError(Exception):

    def __init__(self, feature_line):
        self.feature_line = feature_line

    def __str__(self):
        return "The 'id' property of the feature object on line {fl} is not referenced by a property of any other evidence object.".format(fl=str(self.feature_line))

class DerivedFeatureMismatch(Exception):

    def __init__(self, derived_line, feature_line):
        self.derived_line = derived_line
        self.feature_line = feature_line

    def __str__(self):
        return "The value of the 'derivedFeature' property on line {dl} does not match the value of the 'feature' property on line {fl}.".format(dl=str(self.derived_line),
                                                                                                                                                 fl=str(self.feature_line))

class DerivedFeatureMissing(Exception):

    def __init__(self, alteration_line, feature_line):
        self.alteration_line = alteration_line
        self.feature_line = feature_line
        

    def __str__(self):
        return "The 'derivedFeature' property is missing from the alteration object on line {al}. Its value must match that of the 'feature' property on line {fl}.".format(al=str(self.alteration_line),
                                                                                                                                                                            fl=str(self.feature_line))

class RegionFeatureMissing(Exception):

    def __init__(self, derived_line, region_line):
        self.derived_line = derived_line
        self.region_line = region_line
        
    def __str__(self):
        return "The 'feature' property is missing from the region object on line {rl}. Its value must match that of the 'derivedFeature' property on line {dl}.".format(rl=str(self.region_line),
                                                                                                                                                                        dl=str(self.derived_line))

class DetectionValidator():
    def __init__(self, detection_schema):
        self.detection_schema = detection_schema

    @classmethod
    def __validate_evidence_host_taxa(cls, evidence, evidence_path, detection_map, host_taxa=None, host_path=None):
        if host_taxa and host_path:
            root_host_taxa = host_taxa

            root_host_path = host_path

            if 'hostTaxa' in evidence:
                sub_host_taxa = {host_taxon for host_taxon in evidence['hostTaxa']}

                sub_host_path = evidence_path + ['hostTaxa']

                if len(sub_host_taxa.difference(host_taxa)) > 0:
                    host_position = detection_map['/'.join(host_path)]

                    sub_host_position = detection_map['/'.join(sub_host_path)]

                    raise InconsistentHostTaxa(host_position.key_start.line + 1,
                                               sub_host_position.key_start.line + 1)
        elif 'hostTaxa' in evidence:
            root_host_taxa = {host_taxon for host_taxon in evidence['hostTaxa']}
    
            root_host_path = evidence_path + ['hostTaxa']
        else:
            root_host_taxa = None
    
            root_host_path = None

        return (root_host_taxa, root_host_path)

    @classmethod
    def __validate_feature_host_taxa(cls, feature, feature_path, detection_map, feature_dict, feature_map,
                                     host_taxa=None, host_path=None):
        (root_host_taxa, root_host_path) = cls.__validate_evidence_host_taxa(feature, feature_path, detection_map,
                                                                             host_taxa, host_path)

        if 'subFeatures' in feature:
            for k in range(0, len(feature['subFeatures'])):
                sub_feature_ID = feature['subFeatures'][k]

                if sub_feature_ID in feature_dict:
                    sub_feature = feature_dict[sub_feature_ID]

                    sub_feature_path = ['', 'features', str(feature_map[sub_feature_ID])]

                    cls.__validate_feature_host_taxa(sub_feature, sub_feature_path, detection_map, feature_dict,
                                                     feature_map, host_taxa, host_path)

    @classmethod
    def __validate_alteration_host_taxa(cls, alteration, alteration_path, detection_map, feature_dict, feature_map,
                                        host_taxa=None, host_path=None):
        (root_host_taxa, root_host_path) = cls.__validate_evidence_host_taxa(alteration, alteration_path,
                                                                             detection_map, host_taxa, host_path)

        if 'targetRegion' in alteration:
            target_region = alteration['targetRegion']

            target_region_path = alteration_path + ['targetRegion']

            cls.__validate_region_host_taxa(target_region, target_region_path, detection_map, feature_dict,
                                            feature_map, root_host_taxa, root_host_path)

        if 'targetFeature' in alteration and alteration['targetFeature'] in feature_dict:
            target_feature = feature_dict[alteration['targetFeature']]

            target_feature_path = ['', 'features', str(feature_map[alteration['targetFeature']])]

            cls.__validate_feature_host_taxa(target_feature, target_feature_path, detection_map, feature_dict,
                                             feature_map, root_host_taxa, root_host_path)


    @classmethod
    def __validate_region_host_taxa(cls, region, region_path, detection_map, feature_dict, feature_map, host_taxa=None,
                                    host_path=None):
        (root_host_taxa, root_host_path) = cls.__validate_evidence_host_taxa(region, region_path, detection_map,
                                                                             host_taxa, host_path)

        if 'feature' in region and region['feature'] in feature_dict:
            feature = feature_dict[region['feature']]

            feature_path = region_path + ['feature']

            cls.__validate_feature_host_taxa(feature, feature_path, detection_map, feature_dict, feature_map,
                                             root_host_taxa, root_host_path)

        if 'sourceAlteration' in region:
            alteration = region['sourceAlteration']

            alteration_path = region_path + ['sourceAlteration']

            cls.__validate_alteration_host_taxa(alteration, alteration_path, detection_map, feature_dict, feature_map,
                                                root_host_taxa, root_host_path)

    @classmethod
    def __validate_sequencing_datum_host_taxa(cls, seq_datum, seq_datum_path, detection_map, feature_dict, feature_map,
                                              host_taxa=None, host_path=None):
        (root_host_taxa, root_host_path) = cls.__validate_evidence_host_taxa(seq_datum, seq_datum_path, detection_map,
                                                                             host_taxa, host_path)

        if 'regions' in seq_datum:
            for i in range(0, len(seq_datum['regions'])):
                region = seq_datum['regions'][i]

                region_path = seq_datum_path + ['regions', str(i)]

                cls.__validate_region_host_taxa(region, region_path, detection_map, feature_dict, feature_map,
                                                root_host_taxa, root_host_path)

    @classmethod
    def __validate_detection_host_taxa(cls, detection_data, detection_map, feature_dict, feature_map):
        if 'detections' in detection_data:
            for i in range(0, len(detection_data['detections'])):
                detection = detection_data['detections'][i]

                detection_path = ['', 'detections', str(i)]

                (host_taxa, host_path) = cls.__validate_evidence_host_taxa(detection, detection_path,
                                                                           detection_map)

                if 'reads' in detection:
                    for j in range(0, len(detection['reads'])):
                        read = detection['reads'][j]

                        read_path = detection_path + ['reads', str(j)]

                        cls.__validate_sequencing_datum_host_taxa(read, read_path, detection_map, feature_dict,
                                                                  feature_map, host_taxa, host_path)

                if 'contigs' in detection:
                    for j in range(0, len(detection['contigs'])):
                        contig = detection['contigs'][j]

                        contig_path = detection_path + ['contigs', str(j)]

                        cls.__validate_sequencing_datum_host_taxa(contig, contig_path, detection_map, feature_dict,
                                                                  feature_map, host_taxa, host_path)

    @classmethod
    def __validate_region_coordinates(cls, region, parent_key, parent_indices, detection_map,
                                      alteration_level=0):
        region_path = ['',
                       'detections',
                       str(parent_indices[0]),
                       parent_key,
                       str(parent_indices[1]),
                       'regions',
                       str(parent_indices[2])]

        n = 0
        while (n < alteration_level):
            region_path = region_path + ['sourceAlteration', 'targetRegion']

            n = n + 1

        if region['start'] > region['end']:
            start_position = detection_map['/'.join(region_path + ['start'])]
            end_position = detection_map['/'.join(region_path + ['end'])]

            raise CoordinatesMisordered('start', 'end', start_position.key_start.line + 1,
                                        end_position.key_start.line + 1)

        if (('featureStart' in region or 'featureEnd' in region or 'featureOrientation' in region)
                and not 'feature' in region):
            region_position = detection_map['/'.join(region_path)]

            coord_keys = []

            if 'featureStart' in region:
                coord_keys.append('featureStart')

            if 'featureEnd' in region:
                coord_keys.append('featureEnd')

            if 'featureOrientation' in region:
                coord_keys.append('featureOrientation')

            raise CoordinatesMissingSubject('region', 'feature', region_position.value_start.line + 1, coord_keys)

        if 'featureStart' in region and 'featureEnd' in region:
            feature_start_position = detection_map['/'.join(region_path + ['featureStart'])]
            feature_end_position = detection_map['/'.join(region_path + ['featureEnd'])]

            if region['featureStart'] > region['featureEnd']:
                raise CoordinatesMisordered('featureStart', 'featureEnd',
                                            feature_start_position.key_start.line + 1,
                                            feature_end_position.key_start.line + 1)

            # if 'feature' in region and region['feature'] in feature_dict:
            #     feature_entry = feature_dict[region['feature']]

            #     if region['featureEnd'] - region['featureStart'] + 1 != feature_entry['seqLength']:
            #         raise CoordinatesLengthMismatch('featureStart', 'featureEnd', feature_start_position.key_start.line + 1,
            #                                         feature_end_position.key_start.line + 1, feature_entry['seqLine'])

        if 'sourceAlteration' in region:
            alteration = region['sourceAlteration']

            if 'targetRegion' in alteration:
                target_region = alteration['targetRegion']

                cls.__validate_region_coordinates(target_region, parent_key, parent_indices, detection_map,
                                                  alteration_level + 1)

    @classmethod
    def __validate_detection_coordinates(cls, detection_data, detection_map):
        if 'features' in detection_data:
            for i in range(0, len(detection_data['features'])):
                feature = detection_data['features'][i]

                feature_path = ['', 'features', str(i)]

                if (('sourceStart' in feature or 'sourceEnd' in feature or 'sourceOrientation' in feature)
                        and not 'source' in feature):
                    feature_position = detection_map['/'.join(feature_path)]

                    coord_keys = []

                    if 'sourceStart' in feature:
                        coord_keys.append('sourceStart')

                    if 'sourceEnd' in feature:
                        coord_keys.append('sourceEnd')

                    if 'sourceOrientation' in feature:
                        coord_keys.append('sourceOrientation')

                    raise CoordinatesMissingSubject('feature', 'source', feature_position.value_start.line + 1,
                                                    coord_keys)

                if ('sourceStart' in feature and 'sourceEnd' in feature
                        and feature['sourceStart'] > feature['sourceEnd']):
                    source_start_position = detection_map['/'.join(feature_path + ['sourceStart'])]

                    source_end_position = detection_map['/'.join(feature_path + ['sourceEnd'])]

                    raise CoordinatesMisordered('sourceStart', 'sourceEnd',
                                                source_start_position.key_start.line + 1,
                                                source_end_position.key_start.line + 1)

                # feature_dict[feature['id']] = {}

                # if 'sequence' in feature:
                    # seq_position = detection_map['/'.join(feature_path + ['sequence'])]

                    # feature_dict[feature['id']]['seqLength'] = len(feature['sequence'])
                    # feature_dict[feature['id']]['seqLine'] = seq_position.key_start.line + 1

        if 'detections' in detection_data:
            for i in range(0, len(detection_data['detections'])):
                detection = detection_data['detections'][i]

                if 'reads' in detection:
                    for j in range(0, len(detection['reads'])):
                        read = detection['reads'][j]

                        if 'regions' in read:
                            for k in range(0, len(read['regions'])):
                                region = read['regions'][k]

                                cls.__validate_region_coordinates(region, 'reads', (i, j, k), detection_map)

                if 'contigs' in detection:
                    for j in range(0, len(detection['contigs'])):
                        contig = detection['contigs'][j]

                        if 'regions' in contig:
                            for k in range(0, len(contig['regions'])):
                                region = contig['regions'][k]

                                cls.__validate_region_coordinates(region, 'contigs', (i, j, k), detection_map)

    @classmethod
    def __validate_region_features(cls, region, parent_key, parent_indices, detection_map, feature_IDs):
        feature_refs = set()

        if 'feature' in region:
            feature_refs.add(region['feature'])

            if region['feature'] not in feature_IDs:
                feature_position = detection_map['/'.join(['',
                                                           'detections',
                                                           str(parent_indices[0]),
                                                           parent_key,
                                                           str(parent_indices[1]),
                                                           'regions',
                                                           str(parent_indices[2]),
                                                           'feature'])]

                raise FeatureReferenceError(region['feature'],
                                            'feature',
                                            feature_position.key_start.line + 1)

        if 'sourceAlteration' in region:
            alteration = region['sourceAlteration']

            if 'targetFeature' in alteration:
                feature_refs.add(alteration['targetFeature'])

                if alteration['targetFeature'] not in feature_IDs:
                    target_position = detection_map['/'.join(['',
                                                              'detections',
                                                              str(parent_indices[0]),
                                                              parent_key,
                                                              str(parent_indices[1]),
                                                              'regions',
                                                              str(parent_indices[2]),
                                                              'sourceAlteration',
                                                              'targetFeature'])]

                    raise FeatureReferenceError(alteration['targetFeature'],
                                                'targetFeature',
                                                target_position.key_start.line + 1)

            alteration_position = detection_map['/'.join(['',
                                                          'detections',
                                                          str(parent_indices[0]),
                                                          parent_key,
                                                          str(parent_indices[1]),
                                                          'regions',
                                                          str(parent_indices[2]),
                                                          'sourceAlteration'])]

            region_position = detection_map['/'.join(['',
                                                      'detections',
                                                      str(parent_indices[0]),
                                                      parent_key,
                                                      str(parent_indices[1]),
                                                      'regions',
                                                      str(parent_indices[2])])]

            if 'derivedFeature' in alteration:
                feature_refs.add(alteration['derivedFeature'])

                derived_position = detection_map['/'.join(['',
                                                           'detections',
                                                           str(parent_indices[0]),
                                                           parent_key,
                                                           str(parent_indices[1]),
                                                           'regions',
                                                           str(parent_indices[2]),
                                                           'sourceAlteration',
                                                           'derivedFeature'])]

                if 'feature' not in region:
                    raise RegionFeatureMissing(derived_position.key_start.line + 1,
                                               region_position.value_start.line + 1)
                elif alteration['derivedFeature'] != region['feature']:
                    feature_position = detection_map['/'.join(['',
                                                               'detections',
                                                               str(parent_indices[0]),
                                                               parent_key,
                                                               str(parent_indices[1]),
                                                               'regions',
                                                               str(parent_indices[2]),
                                                               'feature'])]

                    raise DerivedFeatureMismatch(derived_position.key_start.line + 1,
                                                 feature_position.key_start.line + 1)
            elif 'feature' in region:
                feature_position = detection_map['/'.join(['',
                                                           'detections',
                                                           str(parent_indices[0]),
                                                           parent_key,
                                                           str(parent_indices[1]),
                                                           'regions',
                                                           str(parent_indices[2]),
                                                           'feature'])]

                raise DerivedFeatureMissing(alteration_position.key_start.line + 1,
                                            feature_position.key_start.line + 1)

        return feature_refs

    @classmethod
    def __validate_detection_features(cls, detection_data, detection_map, feature_IDs, feature_map):
        feature_refs = set()

        if 'features' in detection_data:
            for i in range(0, len(detection_data['features'])):
                feature = detection_data['features'][i]

                if 'subFeatures' in feature:
                    for sub_feature_ID in feature['subFeatures']:
                        feature_refs.add(sub_feature_ID)

                        if sub_feature_ID not in feature_IDs:
                            sub_feature_position = detection_map['/'.join(['',
                                                                           'features',
                                                                           str(i),
                                                                           'subFeatures'])]

                            raise FeatureReferenceError(sub_feature_ID,
                                                        'subFeatures',
                                                        sub_feature_position.key_start.line + 1)

        if 'detections' in detection_data:
            for i in range(0, len(detection_data['detections'])):
                detection = detection_data['detections'][i]

                if 'reads' in detection:
                    for j in range(0, len(detection['reads'])):
                        read = detection['reads'][j]

                        if 'regions' in read:
                            for k in range(0, len(read['regions'])):
                                region = read['regions'][k]

                                feature_refs.update(cls.__validate_region_features(region, 'reads', (i, j, k),
                                                                                   detection_map, feature_IDs))

                if 'contigs' in detection:
                    for j in range(0, len(detection['contigs'])):
                        contig = detection['contigs'][j]

                        if 'regions' in contig:
                            for k in range(0, len(contig['regions'])):
                                region = contig['regions'][k]

                                feature_refs.update(cls.__validate_region_features(region, 'contigs', (i, j, k),
                                                                                   detection_map, feature_IDs))

                if 'alterations' in detection:
                    for alteration in detection['alterations']:
                        if 'derivedFeature' in alteration:
                            feature_refs.add(alteration['derivedFeature'])

                        if 'targetFeature' in alteration:
                            feature_refs.add(alteration['targetFeature'])

        for feature_ID in feature_IDs:
            if feature_ID not in feature_refs:
                feature_path = ['', 'features', str(feature_map[feature_ID])]

                feature_position = detection_map['/'.join(feature_path)]

                raise FeatureIDError(feature_position.value_start.line + 1)

    def preprocess_and_validate(self, detection_file):
        with open(detection_file) as detection_handle:
            detection_data = json.load(detection_handle)

        with open(detection_file) as detection_handle:
            detection_map = calculate(detection_handle.read())

        self.validate(detection_data, detection_map)
            
    def validate(self, detection_data, detection_map):
        validate(instance=detection_data, schema=self.detection_schema)

        if 'features' in detection_data:
            feature_dict = {feature['id'] : feature for feature in detection_data['features']}
            feature_map = {detection_data['features'][i]['id'] : i for i in range(0, len(detection_data['features']))}
        else:
            feature_dict = {}
            feature_map = {}

        self.__validate_detection_host_taxa(detection_data, detection_map, feature_dict, feature_map)

        self.__validate_detection_features(detection_data, detection_map, feature_dict.keys(), feature_map)

        self.__validate_detection_coordinates(detection_data, detection_map)

def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--schema_file')
    parser.add_argument('-d', '--detection_files', nargs='*', default=[])
    
    args = parser.parse_args(args)

    with open(args.schema_file) as schema_handle:
        detection_schema = json.load(schema_handle)

    detection_validator = DetectionValidator(detection_schema)

    detection_files = []
    for detection_file in args.detection_files:
        if os.path.isdir(detection_file):
            detection_files.extend([os.path.join(detection_file, df) for df in os.listdir(detection_file) if
                                    os.path.isfile(os.path.join(detection_file, df)) and df.endswith('.json')])
        else:
            detection_files.append(detection_file)

    for detection_file in detection_files:
        print('Validating ' + detection_file + '...')

        with open(detection_file) as detection_handle:
            detection_data = json.load(detection_handle)

        with open(detection_file) as detection_handle:
            detection_map = calculate(detection_handle.read())

        detection_validator.validate(detection_data, detection_map)
        
        print('Determined ' + detection_file + ' is valid.')


if __name__ == '__main__':
    main()