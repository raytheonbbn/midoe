#!/usr/bin/env python3
import copy
import logging
import os
import sys
import argparse
import json
import pandas as pd
from schema.dataschema import Detection, Investigation, dumpJSON, jsonobjproc
from evidenceintegration.DetectionsToVis import DetectionsToVisualization
from evidenceintegration.EvidenceIntegration import EvidenceIntegration, RegionGroup, AlterationGroup
from evidenceintegration.PreprocessInvestigations import PreprocessInvestigations
from evidenceintegration.utils import EvidenceIntegrationUtils

class DetectionsToCSV:

    NOT_APPLICABLE = "Not applicable"
    TODO = "TODO"
    SAMPLE_ID = "Sample_ID"
    NUM_ASSAYS = "Number_of_Technical_and/or_Biological_Replicates_of_Assay_Performed"
    CONCENTRATION = "Concentration_and_Volume_Per_Assay"
    REF_GEN_USED = "Reference_Genome_or_Assembly_Used (Y/N)"
    REF_NAME_ACC = "Reference_Name_or_Accession"
    DECISION_CRITERIA = "Decision_Criteria"
    DNA_CONCENTRATION = "Starting_DNA_Concentration"
    DNA_CONCENTRATION_UNITS = "Starting_DNA_Concentration_Units"
    DATE_RUN = "Date_Run"
    TESTER = "Tester_Name"
    SW_VER = "Software_Version"
    HW_ENV = "Hardware_Environment"
    COMP_RES = "Computational_Resources"
    COMP_TIME = "Computational_Analysis_Time"
    ASSAY_TIME = "Experimental_Assay_Time"
    SCORE = "Score_or_Confidence_Measure"
    ENG_DETECTED = "Engineering_Detected"
    ENG_NAME = "What_was_Detected"
    NATIVE = "Native_to_Host"
    EVIDENCE = "Evidence_of_Engineering"
    HOST_SPECIES = "Host_Species"
    ENG_DNA_SOURCE = "Source_of_Engineered_DNA"
    COORDINATES = "Base_Pair_Coordinates_or_Gene_Context_in_the_T_&_E_sample"
    SIZE = "Size_of_Signature"
    C_OR_P = "Chromosome_or_Plasmid"
    CHROMOSOME = "Which_Chromosome"
    PART_CLASS = "Part_Class"
    DETECTION_MODULE = "Detection_Module"
    NOTES = "Notes"
    SIGNATURE_ID = "Signature_ID"
    SIGNATURE_GROUP = "Signature_Group"
    PARENT_SIGNATURE = "Parent_Signature"
    READ_COUNT = "Read_Count"
    SEQUENCE = "Signature_Sequence"    

    MIN_CSV_COLUMNS = [
        SAMPLE_ID,
        REF_NAME_ACC,
        SCORE,
        ENG_DETECTED,
        ENG_NAME,
        NATIVE,
        EVIDENCE,
        HOST_SPECIES,
        COORDINATES,
        SIZE,
        PART_CLASS,
        DETECTION_MODULE,
        SIGNATURE_ID,
        SIGNATURE_GROUP,
        READ_COUNT,
        SEQUENCE
    ]

    TE_CSV_COLUMNS = [
        SAMPLE_ID,
        NUM_ASSAYS,
        CONCENTRATION,
        REF_GEN_USED,
        REF_NAME_ACC,
        DECISION_CRITERIA,
        DNA_CONCENTRATION,
        DNA_CONCENTRATION_UNITS,
        DATE_RUN,
        TESTER,
        SW_VER,
        HW_ENV,
        COMP_RES,
        COMP_TIME,
        ASSAY_TIME,
        SCORE,
        ENG_DETECTED,
        ENG_NAME,
        NATIVE,
        EVIDENCE,
        HOST_SPECIES,
        ENG_DNA_SOURCE,
        COORDINATES,
        SIZE,
        C_OR_P,
        CHROMOSOME,
        PART_CLASS,
        DETECTION_MODULE,
        NOTES,
        SIGNATURE_ID,
        SIGNATURE_GROUP,
        PARENT_SIGNATURE,
        READ_COUNT,
        SEQUENCE
    ]


    @staticmethod
    def natural_te_row_for_sample(sample_id):
        row_dict = {}
        row_dict[DetectionsToCSV.SAMPLE_ID] = sample_id
        row_dict[DetectionsToCSV.NUM_ASSAYS] = 1
        row_dict[DetectionsToCSV.CONCENTRATION] = DetectionsToCSV.NOT_APPLICABLE
        row_dict[DetectionsToCSV.REF_GEN_USED] = "no"
        row_dict[DetectionsToCSV.REF_NAME_ACC] = DetectionsToCSV.NOT_APPLICABLE
        row_dict[DetectionsToCSV.DECISION_CRITERIA] = DetectionsToCSV.TODO
        row_dict[DetectionsToCSV.DNA_CONCENTRATION] = DetectionsToCSV.NOT_APPLICABLE
        row_dict[DetectionsToCSV.DNA_CONCENTRATION_UNITS] = DetectionsToCSV.NOT_APPLICABLE
        row_dict[DetectionsToCSV.DATE_RUN] = ""
        row_dict[DetectionsToCSV.TESTER] = DetectionsToCSV.TODO
        row_dict[DetectionsToCSV.SW_VER] = DetectionsToCSV.TODO
        row_dict[DetectionsToCSV.HW_ENV] = DetectionsToCSV.TODO
        row_dict[DetectionsToCSV.COMP_RES] = DetectionsToCSV.TODO
        row_dict[DetectionsToCSV.COMP_TIME] = DetectionsToCSV.TODO
        row_dict[DetectionsToCSV.ASSAY_TIME] = DetectionsToCSV.NOT_APPLICABLE
        row_dict[DetectionsToCSV.SCORE] = "0"
        row_dict[DetectionsToCSV.ENG_DETECTED] = "no"
        row_dict[DetectionsToCSV.ENG_NAME] = DetectionsToCSV.NOT_APPLICABLE
        row_dict[DetectionsToCSV.NATIVE] = DetectionsToCSV.NOT_APPLICABLE
        row_dict[DetectionsToCSV.EVIDENCE] = DetectionsToCSV.NOT_APPLICABLE
        row_dict[DetectionsToCSV.HOST_SPECIES] = DetectionsToCSV.TODO
        row_dict[DetectionsToCSV.ENG_DNA_SOURCE] = DetectionsToCSV.NOT_APPLICABLE
        row_dict[DetectionsToCSV.COORDINATES] = DetectionsToCSV.NOT_APPLICABLE
        row_dict[DetectionsToCSV.SIZE] = DetectionsToCSV.NOT_APPLICABLE
        row_dict[DetectionsToCSV.C_OR_P] = DetectionsToCSV.NOT_APPLICABLE
        row_dict[DetectionsToCSV.CHROMOSOME] = DetectionsToCSV.NOT_APPLICABLE
        row_dict[DetectionsToCSV.PART_CLASS] = DetectionsToCSV.NOT_APPLICABLE
        row_dict[DetectionsToCSV.DETECTION_MODULE] = "GUARDIAN"
        row_dict[DetectionsToCSV.NOTES] = ""

        # -----
        # extra GUARDIAN fields
        # -----
        row_dict[DetectionsToCSV.SIGNATURE_ID] = DetectionsToCSV.NOT_APPLICABLE
        row_dict[DetectionsToCSV.SIGNATURE_GROUP] = DetectionsToCSV.NOT_APPLICABLE
        row_dict[DetectionsToCSV.PARENT_SIGNATURE] = DetectionsToCSV.NOT_APPLICABLE
        row_dict[DetectionsToCSV.SEQUENCE] = DetectionsToCSV.NOT_APPLICABLE


        rows = [row_dict]
        df = pd.DataFrame(rows,
            columns=DetectionsToCSV.MIN_CSV_COLUMNS)
        return df


    @staticmethod
    def investigation_to_te_df(investigation):
        detections_of_interest = []
        df = pd.DataFrame(
            columns=DetectionsToCSV.MIN_CSV_COLUMNS)
        rows = []



        for detection in investigation.get_detections():
            if detection.get_agent() == EvidenceIntegration.GUARDIAN_FEATURE_GROUP_AGENT_TAG:
                continue

            detections_of_interest.append(detection)

        all_features = investigation.get_features()

        # logging.debug("in investigation_to_te_df, all_features: %s", str([f.get_id() for f in all_features]))
        # logging.debug("investigation_to_te_df, investigation provided has %d detections_of_interest", len(detections_of_interest))
        for detection in detections_of_interest:

            # alterations only create one row
            for alteration in detection.get_alterations():
                alteration_rows = DetectionsToCSV.te_rows_for_alteration(alteration, detection, all_features)
                # logging.debug("Appending row for alteration: %s", str(alteration_rows))
                rows = rows + alteration_rows
                # rows.append(row)

            # reads can have multiple regions which point at their own features...
            # 1...n rows
            for read in detection.get_reads():
                # logging.debug("read: %s", str(read))
                read_rows = DetectionsToCSV.te_rows_for_read(read, detection, all_features)
                # logging.debug("Appending rows for read: %s", str(read_rows))
                # rows.append(read_rows)
                rows = rows + read_rows

            # contigs can have multiple regions which point at their own features...
            # 1...n rows
            for contig in detection.get_contigs():
                contig_rows = DetectionsToCSV.te_rows_for_contig(contig, detection, all_features)
                # logging.debug("Appending rows for contig: %s", str(contig_rows))
                # rows.append(contig_rows)
                rows = rows + contig_rows
                



        df = pd.DataFrame(
            rows, columns=DetectionsToCSV.MIN_CSV_COLUMNS)
        return df


    @staticmethod
    def te_rows_for_alteration(alteration, detection, features):
        row_dicts = []

        row_dict = DetectionsToCSV.te_base_row_for_detection(detection)
        # logging.debug("row_dict after te_base_row_for_detection: %s", str(row_dict))

        # row_dict[ENG_NAME] =    # must be set by read, contig, or alteration
        # row_dict[EVIDENCE] =   # must be set by read, contig, or alteration
        # row_dict[DETECTION_MODULE] =  # must be set by read, contig, or alteration


        # overwrite the type, always found in alteration so no try block necessary ATM
        row_dict[DetectionsToCSV.EVIDENCE] = alteration.get_type()

        # try to grab agent
        try:
            if alteration.get_agent() is not None and alteration.get_agent() != "":
                row_dict[DetectionsToCSV.DETECTION_MODULE] = alteration.get_agent()
        except AttributeError:
            # no agent in the alteration... we'll have to get it from feature.
            pass


        # try to grab hostTaxa
        # todo maybe move this into something for all subclasses of Evidence
        try:
            if alteration.get_hostTaxa() is not None and alteration.get_hostTaxa() != "":
                row_dict[DetectionsToCSV.HOST_SPECIES] = alteration.get_hostTaxa()
        except AttributeError:
            pass

        # 2. The "Reference_Genome_or_Assembly_Used_(Y/N)" column should be set to "yes"
        # and "Reference_Name_or_Accession" should contain an accession number 
        # if an Alteration has a "targetAssembly" property.
        try:
            if alteration.get_targetAssembly() is not None:
                row_dict[DetectionsToCSV.REF_GEN_USED] = "yes"
                row_dict[DetectionsToCSV.REF_NAME_ACC] = alteration.get_targetAssembly()
        except AttributeError:
            pass

        # get readCount if it exists (only Targeted Search ATM)..
        if alteration.get_readCount() is not None:
            row_dict[DetectionsToCSV.READ_COUNT] = alteration.get_readCount()


        if alteration.get_derivedFeature() is not None:
            # logging.debug("alteration's derivedFeature: %s", str(alteration.get_derivedFeature()))

            for feature in features:
                if feature.get_id() == alteration.get_derivedFeature():
                    row_dicts = row_dicts + DetectionsToCSV.te_rows_for_feature(feature, row_dict, features)
                    break

        if (alteration.get_derivedFeature() is None) and (alteration.get_targetFeature() is not None):
            # logging.debug("alteration's targetFeature: %s", str(alteration.get_targetFeature()))
            for feature in features:
                if feature.get_id() == alteration.get_targetFeature():
                    row_dicts = row_dicts + DetectionsToCSV.te_rows_for_feature(feature, row_dict, features)
                    break


        return row_dicts


    @staticmethod
    def te_rows_for_contig(contig, detection, features):
        row_dicts = []

        # watch this.. but if the contig doesn't have a region.. we will not produce a T&E row for it.
        regions = []
        try:
            regions = contig.get_regions()
        except AttributeError:
            pass

        for region in regions:
             row_dict = DetectionsToCSV.te_base_row_for_detection(detection)
             row_dicts = row_dicts + DetectionsToCSV.te_rows_for_region(region, row_dict, contig.get_source(), features)
             # row_dicts.append(row_dict)

        return row_dicts


    @staticmethod
    def te_rows_for_read(read, detection, features):
        row_dicts = []

        for region in read.get_regions():
             row_dict = DetectionsToCSV.te_base_row_for_detection(detection)
             row_dicts = row_dicts + DetectionsToCSV.te_rows_for_region(region, row_dict, read.get_source(), features)
             # row_dicts.append(row_dict)

        return row_dicts


    @staticmethod
    def te_rows_for_feature(feature, row_dict_arg, features, parent_feature_id = None):
        row_dicts = []
        row_dict = copy.deepcopy(row_dict_arg)
        # logging.debug("te_rows_for_feature called, feature: %s", str(feature))

        # only set PART_CLASS to Role if it wasn't set by a parent..
        if row_dict[DetectionsToCSV.PART_CLASS] == "":
            try:
                row_dict[DetectionsToCSV.PART_CLASS] = feature.get_role()
            except AttributeError:
                pass


        # get Agent from Feature if it exists
        if feature.get_agent() is not None and feature.get_agent() != "":
            row_dict[DetectionsToCSV.DETECTION_MODULE] = feature.get_agent()


        # "Signature_ID" should be populated with the "id" property of a feature.
        if feature.get_id() is not None and feature.get_id() != "":
            row_dict[DetectionsToCSV.SIGNATURE_ID] = feature.get_id()

        # "Signature_Sequence" should be populated the "sequence" property of the feature.
        if feature.get_sequence() is not None and feature.get_sequence() != "":
            row_dict[DetectionsToCSV.SEQUENCE] = feature.get_sequence()
            row_dict[DetectionsToCSV.SIZE] = len(feature.get_sequence())

        # "Parent_Signature" should be populated if a row represents a feature identified
        # by the "subFeatures" property of another feature (populate using this parent feature's
        # "id" property).
        if parent_feature_id is not None:
            row_dict[DetectionsToCSV.PARENT_SIGNATURE] = parent_feature_id

        identifier = ""
        feature_id = feature.get_id()
        feature_name = feature.get_name()
        feature_source = feature.get_source()


        # If a Feature has both a name and a source, compare them. If either its name
        # or its source is a substring of the other, then identify the feature using its source.
        # Otherwise, identify the feature using the concatenation of “Subsequence of “ and 
        # its name, stripping any commas from the result.
        if (feature_name is not None) and (feature_source is not None):
            if (feature_name in feature_source) or (feature_source in feature_name):
                identifier = feature_source
            else:
                identifier = "Subsequence of " + feature_name
                identifier = "".join(identifier.split(","))

        # If a Feature lacks a name but has a source, then identify the Feature using the concatenation
        # of of “Subsequence of “ and its source, stripping any commas from the result.
        elif (feature_name is None) and (feature_source is not None):
            identifier = "Subsequence of " + feature_source
            identifier = "".join(identifier.split(","))

        # If a Feature lacks a source but has a name, then identify the Feature using its name.
        elif (feature_source is None) and (feature_name is not None):
            identifier = feature_name

        # If a Feature lacks both a source and a name, then identify it using its id.
        elif (feature_source is None) and (feature_name is None):
            identifier = feature_id

        row_dict[DetectionsToCSV.ENG_NAME] = identifier
        row_dicts.append(row_dict)

        # If a Feature has subFeatures, then apply the previous steps to them as well. 
        # For the overall JSON-to-CSV conversion, this should yield one row for the parent 
        # Feature and one row for each child Feature. 
        if (feature.get_subFeatures() is not None) and (len(feature.get_subFeatures()) > 0):
            for feature_in_subfeatures in feature.get_subFeatures():
                for f in features:
                    if f.get_id() == feature_in_subfeatures:
                        row_dicts = row_dicts + DetectionsToCSV.te_rows_for_feature(f, row_dict, features, feature.get_id())
                        # row_dicts = row_dicts + DetectionsToCSV.te_rows_for_feature(feature, row_dict_arg, features)



        return row_dicts


    @staticmethod
    def te_rows_for_region(region, row_dict_arg, source_arg, features):
        row_dicts = []
        row_dict = copy.deepcopy(row_dict_arg)
        # logging.debug("te_rows_for_region called, region: %s", str(region))


        # source, start, and end
        source = source_arg.split("/")[-1]
        start = -1
        end = -1
        try:
            start = region.get_start()
            end = region.get_end()
        except AttributeError:
            pass

        try:
            start = region.get_featureStart()
            end = region.get_featureEnd()
        except AttributeError:
            pass

        if start != -1 and end != -1:
            coordinates = source + ":" + str(start) + "-" + str(end)
        else:
            coordinates = source

        row_dict[DetectionsToCSV.COORDINATES] = coordinates


        # role into Part_Class
        try:
            row_dict[DetectionsToCSV.PART_CLASS] = region.get_role()
        except AttributeError:
            pass


        # get info from the region's feature if it has one (it better...)
        if region.get_feature() is not None:
            for feature in features:
                if feature.get_id() == region.get_feature():
                    row_dicts = row_dicts + DetectionsToCSV.te_rows_for_feature(feature, row_dict, features)


        return row_dicts



    @staticmethod
    def te_base_row_for_detection(detection):

        row_dict = {}
        row_dict[DetectionsToCSV.SAMPLE_ID] = detection.get_sample()
        row_dict[DetectionsToCSV.NUM_ASSAYS] = 1
        row_dict[DetectionsToCSV.CONCENTRATION] = DetectionsToCSV.NOT_APPLICABLE
        row_dict[DetectionsToCSV.REF_GEN_USED] = "no"
        row_dict[DetectionsToCSV.REF_NAME_ACC] = DetectionsToCSV.NOT_APPLICABLE
        row_dict[DetectionsToCSV.DECISION_CRITERIA] = DetectionsToCSV.TODO
        row_dict[DetectionsToCSV.DNA_CONCENTRATION] = DetectionsToCSV.NOT_APPLICABLE
        row_dict[DetectionsToCSV.DNA_CONCENTRATION_UNITS] = DetectionsToCSV.NOT_APPLICABLE
        try:
            row_dict[DetectionsToCSV.DATE_RUN] = detection.get_date()
        except AttributeError:
            row_dict[DetectionsToCSV.DATE_RUN] = ""
        row_dict[DetectionsToCSV.TESTER] = DetectionsToCSV.TODO
        row_dict[DetectionsToCSV.SW_VER] = DetectionsToCSV.TODO
        row_dict[DetectionsToCSV.HW_ENV] = DetectionsToCSV.TODO
        row_dict[DetectionsToCSV.COMP_RES] = DetectionsToCSV.TODO
        row_dict[DetectionsToCSV.COMP_TIME] = DetectionsToCSV.TODO
        row_dict[DetectionsToCSV.ASSAY_TIME] = DetectionsToCSV.NOT_APPLICABLE
        try:
            row_dict[DetectionsToCSV.SCORE] = detection.get_confidence()
        except AttributeError:
            row_dict[DetectionsToCSV.SCORE] = ""
        row_dict[DetectionsToCSV.ENG_DETECTED] = "yes"
        row_dict[DetectionsToCSV.ENG_NAME] = ""   # must be set by read, contig, or alteration
        row_dict[DetectionsToCSV.NATIVE] = "no"
        row_dict[DetectionsToCSV.EVIDENCE] = "insertion"
        row_dict[DetectionsToCSV.HOST_SPECIES] = ""
        row_dict[DetectionsToCSV.ENG_DNA_SOURCE] = ""
        row_dict[DetectionsToCSV.COORDINATES] = ""
        row_dict[DetectionsToCSV.SIZE] = ""
        row_dict[DetectionsToCSV.C_OR_P] = DetectionsToCSV.TODO
        row_dict[DetectionsToCSV.CHROMOSOME] = DetectionsToCSV.TODO
        row_dict[DetectionsToCSV.PART_CLASS] = ""
        row_dict[DetectionsToCSV.DETECTION_MODULE] = "" # must be set by read, contig, or alteration
        row_dict[DetectionsToCSV.NOTES] = ""

        # -----
        # extra GUARDIAN fields
        # -----
        row_dict[DetectionsToCSV.SIGNATURE_ID] = ""
        row_dict[DetectionsToCSV.SIGNATURE_GROUP] = detection.get_sample() + ":" +  DetectionsToCSV.signature_group_str_for_detection(detection)
        row_dict[DetectionsToCSV.PARENT_SIGNATURE] = ""
        row_dict[DetectionsToCSV.READ_COUNT] = ""
        row_dict[DetectionsToCSV.SEQUENCE] = ""

        return row_dict


    @staticmethod
    def signature_group_str_for_detection(detection):
        # For the "Signature_Group" column, we could potentially use any appropriate IDs for
        # the meta-groups that you may form during integration, or we could use the shortest
        # ID of one of the features from a group. For the attached example, I made the feature
        # IDs by prefixing the pre-evidence-integration feature IDs with their detecting agents,
        # but these IDs would actually take whatever form we are currently using to ensure unique
        # feature IDs post-evidence-integration.

        # logging.debug("Detection, singular=%d, agent=%s", int(detection.get_singular()), str(detection.get_agent()))

        # These are from EvidenceIntegration, so all_regions are in one RegionGroup
        # and all_alterations are in one AlterationGroup
        all_regions = EvidenceIntegration.all_regions_for_detection(detection)
        all_alterations = EvidenceIntegration.all_alterations_for_detection(
            detection)

        alteration_str = None
        # it has an AlterationGroup
        if len(all_alterations):
            alteration_group = AlterationGroup(all_alterations[0])
            for alteration in all_alterations[1:]:
                alteration_group.add_alteration(alteration)

            try:
                alteration_str = alteration_group.descriptive_name_for_csv()
            except AttributeError as ae:
                logging.warn("Using str(AlterationGroup) due to AttributeError: %s", str(ae))
                alteration_str = str(alteration_group)

        region_str = None
        # it has a RegionGroup
        if len(all_regions):
            region_group = RegionGroup(all_regions[0])
            for region in all_regions[1:]:
                region_group.add_region(region)


            try:
                region_str = region_group.descriptive_name_for_csv()
            except AttributeError as ae:
                logging.warn(
                    "Using str(RegionGroup) due to AttributeError: %s", str(ae))
                region_str = str(region_group)


        if alteration_str is not None and region_str is not None:
            signature_group_str = "Composite detection with an AlterationGroup and RegionGroup. AlterationGroup: "
            signature_group_str += alteration_str + ". RegionGroup: "
            signature_group_str += region_str
        elif alteration_str is not None:
            signature_group_str = "AlterationGroup: "
            signature_group_str += alteration_str
        else:
            signature_group_str = "RegionGroup: "
            signature_group_str += region_str

        # logging.debug("signature_group_str: %s", str(signature_group_str))
        # logging.debug("alteration_str: %s", str(alteration_str))
        # logging.debug("region_str: %s", str(region_str))
        return signature_group_str


def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser()
    parser.add_argument('-j', '--json_files', nargs='*', default=[])
    parser.add_argument('-d', '--json_dirs', nargs='*', default=[])
    parser.add_argument('-o', '--out_dir', nargs=1, default=None)
    parser.add_argument('-l', '--integration_log', nargs='?', default='')
    parser.add_argument('-i', '--integrate', action='store_true', default=False)

    args = parser.parse_args(args)
    json_files = args.json_files
    json_dirs = args.json_dirs
    out_dir = args.out_dir
    integrate = args.integrate

    if len(args.integration_log) > 0:
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)

        logging.basicConfig(level=logging.INFO, filename=args.integration_log, filemode='w',
                            format='%(levelname)s : %(message)s')
    else:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(levelname)s : %(message)s')

    if ((len(json_files) == 0 or json_files is None) and (len(json_dirs) == 0 or json_dirs is None)):
        logging.error(
            "Must supply either a directory in json_dirs, or a file in json_files.")
        exit()

    logging.debug("json_files in: %s", str(json_files))
    if len(json_files) == 1:
        parts = json_files[0].split(",")
        tmp = []
        for part in parts:
            tmp.append(part.strip())
        json_files = tmp

    if (len(json_dirs) == 1):
        parts = json_dirs[0].split(",")
        tmp = []
        for part in parts:
            tmp.append(part.strip())
        for json_dir in tmp:
            json_files = json_files + \
                EvidenceIntegrationUtils.get_json_files_from_dir(
                    json_dir)

    if (len(out_dir) == 1):
        out_dir = out_dir[0]
        os.makedirs(out_dir, exist_ok=True)
        # if not os.path.exists(out_dir):
        #     os.makedirs(out_dir)
    else:
        logging.error("supply output dir for json and csv!")
        exit()

    all_investigations = []
    for json_file in json_files:
        logging.debug(
            "Creating Investigation object for json_file: %s", str(json_file))
        investigation = Investigation()
        with open(json_file) as json_handle:
            json_dict = json.load(json_handle)
            investigation.read_from_json(json_dict)
            logging.debug("investigation: %s", str(investigation))
            all_investigations.append(investigation)



    if integrate:
        base_detections = DetectionsToVisualization.json_files_to_detections(
            json_files)

        frames = []
        processed_investigations = PreprocessInvestigations.process_investigations(
            all_investigations)

        sample_to_guardian_investigation = EvidenceIntegration.integrate(
            processed_investigations, json_out_dir=out_dir)

        logging.debug("all samples: %s", str(sample_to_guardian_investigation.keys()))

        for sample in sample_to_guardian_investigation:
            guardian_investigation = sample_to_guardian_investigation[sample]
            te_df = DetectionsToCSV.investigation_to_te_df(guardian_investigation)
            frames.append(te_df)

        non_eng_sample_ids = []
        for detection in base_detections:
            # skip GUARDIAN detections and GUARDIAN_FEATURE_GROUPS
            try:
                if detection.get_singular() == True:
                    continue
            except AttributeError:
                # doesn't have singular set, go ahead and let it through
                pass
            # if the detection did not result in a guardian investigation, there is no engineering found
            if detection.get_sample() not in sample_to_guardian_investigation:
                te_df = DetectionsToCSV.natural_te_row_for_sample(detection.get_sample())
                non_eng_sample_ids.append(detection.get_sample())
                frames.append(te_df)

        logging.debug("samples w/o engineering found: %s", str(non_eng_sample_ids))

        whole_df = pd.concat(frames, ignore_index=True)
        whole_df.to_csv(out_dir + "/guardian-out.csv", index=False)

    else:
        # just doing 1 for now...
        df = DetectionsToCSV.investigation_to_te_df(all_investigations[0])

        logging.debug("df: %s", str(df))
        df.to_csv(out_dir + "/guardian-out.csv", index=False)

    print('Finished')


if __name__ == '__main__':
    main()

