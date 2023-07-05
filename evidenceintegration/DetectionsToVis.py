#!/usr/bin/env python3
import logging
import os
import json
import pandas as pd
from schema.dataschema import Detection, Investigation, dumpJSON, jsonobjproc
from evidenceintegration.EvidenceIntegration import EvidenceIntegration, RegionGroup, AlterationGroup
import flask


NATURAL = "http://guardian.bbn.technology#natural"
INDETERMINATE = "http://guardian.bbn.technology#indeterminate_origin"
ENGINEERED = "http://guardian.bbn.technology#engineered"
ENGINEERING_PRETTY = {NATURAL: "No",
                      INDETERMINATE: "Indeterminate",
                      ENGINEERED: "Yes"}

TAXON_URI_TO_NAME = {'http://purl.obolibrary.org/obo/NCBITaxon_4932': 'Saccharomyces cerevisiae',
                     'http://purl.obolibrary.org/obo/NCBITaxon_559292': 'Saccharomyces cerevisiae strain S288C',
                     'http://purl.obolibrary.org/obo/NCBITaxon_889517': 'Saccharomyces cerevisiae strain CEN.PK113-7D',
                     'http://purl.obolibrary.org/obo/NCBITaxon_1247190': 'Saccharomyces cerevisiae strain BY4741',
                     'http://purl.obolibrary.org/obo/NCBITaxon_580240': 'Saccharomyces cerevisiae strain W303',
                     'http://purl.obolibrary.org/obo/NCBITaxon_303': 'Pseudomonas putida',
                     'http://purl.obolibrary.org/obo/NCBITaxon_160488': 'Pseudomonas putida strain KT2440',
                     'http://purl.obolibrary.org/obo/NCBITaxon_1423': 'Bacillus subtilis',
                     'http://purl.obolibrary.org/obo/NCBITaxon_224308': 'Bacillus subtilis subspecies subtilis strain 168',
                     'http://purl.obolibrary.org/obo/NCBITaxon_691': 'Vibrio natriegens',
                     'http://purl.obolibrary.org/obo/NCBITaxon_1219067': 'Vibrio natriegens strain ATCC 14048',
                     'http://purl.obolibrary.org/obo/NCBITaxon_562': 'Escherichia coli',
                     'http://purl.obolibrary.org/obo/NCBITaxon_511145': 'Escherichia coli strain K-12 substrain MG1655',
                     'http://purl.obolibrary.org/obo/NCBITaxon_3702': 'Arabidopsis thaliana ecotype Columbia',
                     'http://purl.obolibrary.org/obo/NCBITaxon_287': 'Pseudomonas aeruginosa',
                     'http://purl.obolibrary.org/obo/NCBITaxon_208964': 'Pseudomonas aeruginosa strain PAO1',
                     'http://purl.obolibrary.org/obo/NCBITaxon_546': 'Citrobacter freundii',
                     'http://purl.obolibrary.org/obo/NCBITaxon_2077147': 'Citrobacter freundii complex species CFNIH3',
                     'http://purl.obolibrary.org/obo/NCBITaxon_2697049': 'SARS-CoV-2',
                     'http://purl.obolibrary.org/obo/NCBITaxon_11646': 'Lentivirus',
                     'http://purl.obolibrary.org/obo/NCBITaxon_11676': 'HIV-1',
                     'http://purl.obolibrary.org/obo/NCBITaxon_11320': 'Influenza A virus',
                     'http://purl.obolibrary.org/obo/NCBITaxon_2681598': 'Escherichia phage T4',
                     'http://purl.obolibrary.org/obo/NCBITaxon_10760': 'Escherichia phage T7',
                     'http://purl.obolibrary.org/obo/NCBITaxon_1977402': 'Escherichia virus M13'}


class DetectionsToVisualization:

    # you can remove from this list freely,
    # but if you change these, make sure that "if statements" in the
    # for column in df.columns: loop are corrected
    SAMPLE_OVERVIEW_COLUMNS = ["Batch", "SAMPLE ID", "Engineering Detected", "Consensus",
                               "Confidence", "Host Species", "Detection Module", "Date Run"]

    @staticmethod
    def json_files_to_investigations(json_files):
        # logging.debug(
        #     "json_files_to_overview_table called with json_files: %s", str(json_files))

        all_investigations = []
        for json_file in json_files:
            investigation = Investigation()
            with open(json_file) as json_handle:
                json_dict = json.load(json_handle)
                investigation.read_from_json(json_dict)
                all_investigations.append(investigation)

        return all_investigations

    @staticmethod
    def json_files_to_detections(json_files):
        all_investigations = DetectionsToVisualization.json_files_to_investigations(
            json_files)
        all_detections = []

        for investigation in all_investigations:
            all_detections = all_detections + investigation.get_detections()

        return all_detections

    @staticmethod
    def integration_result_to_detections(sample_to_guardian_investigation):
        all_detections = []

        for sample in sample_to_guardian_investigation:
            try:
                all_detections = all_detections + \
                    sample_to_guardian_investigation[sample].get_detections()
            except AttributeError:
                continue

        return all_detections

    @staticmethod
    def detections_to_sample_overview_df(all_detections):
        detections_of_interest = []
        df = pd.DataFrame(
            columns=DetectionsToVisualization.SAMPLE_OVERVIEW_COLUMNS)
        rows = []

        eng2 = "<img src='" + \
            flask.url_for('static', filename='2eng.png') + "' height='20' />"
        eng1 = "<img src='" + \
            flask.url_for('static', filename='1eng.png') + "' height='20' />"
        eng0 = "<img src='" + \
            flask.url_for('static', filename='0eng.png') + "' height='20' />"
        qmrk = "<img src='" + \
            flask.url_for('static', filename='indeterminate.png') + \
            "' height='20' />"

        for detection in all_detections:
            if detection.get_agent() == EvidenceIntegration.GUARDIAN_FEATURE_GROUP_AGENT_TAG:
                continue

            detections_of_interest.append(detection)

        for detection in detections_of_interest:
            detection_dict = {}
            for column in df.columns:
                if column == "Batch":
                    detection_dict[column] = detection.get_batch()
                if column == "SAMPLE ID":
                    # [Link text Here](https://link-url-here.org)
                    detection_dict[column] = "[" + detection.get_sample() + \
                        "](/analysis?sample_id=" + detection.get_sample() + ")"

                if column == "Engineering Detected":
                    detection_dict[column] = ENGINEERING_PRETTY[
                        detection.get_engineered()
                    ]
                if column == "Consensus":

                    # redo this once scores are available and consistent
                    if detection.get_agent() == EvidenceIntegration.GUARDIAN_AGENT_TAG:
                        detection_dict[column] = eng2
                    else:
                        detection_dict[column] = eng1

                if column == "Confidence":
                    try:
                        detection_dict[column] = detection.get_confidence()
                    except AttributeError:
                        logging.debug(
                            "Detection does not have a confidence value, setting it to 0.")
                        detection_dict[column] = 0
                if column == "Host Species":
                    try:
                        # TODO BBAS: this [0] on the following line is sketchy..
                        try:
                            detection_dict[column] = TAXON_URI_TO_NAME[detection.get_hostTaxa()[
                                0]]
                        except KeyError:
                            detection_dict[column] = detection.get_hostTaxa()
                    except AttributeError:
                        logging.debug(
                            "Detection does not have a hostTaxa, setting it to 'Not Found in source JSON'")
                        detection_dict[column] = "Not Found in source JSON"
                if column == "Detection Module":
                    detection_dict[column] = detection.get_agent()
                try:
                    if column == "Date Run":
                        detection_dict[column] = detection.get_date()
                except AttributeError:
                    detection_dict[column] = None

            # logging.debug("detection_dict: %s", str(detection_dict))
            rows.append(detection_dict)

        df = pd.DataFrame(
            rows, columns=DetectionsToVisualization.SAMPLE_OVERVIEW_COLUMNS)
        return df

    @staticmethod
    def name_and_scores_for_composite_detection(detection, feature_groups=None):
        # logging.debug(
        #     "Making name_and_scores_for_composite_detection: %s", str(detection))

        # These are from EvidenceIntegration, so all_regions are in one RegionGroup
        # and all_alterations are in one AlterationGroup
        all_regions = EvidenceIntegration.all_regions_for_detection(detection)
        all_alterations = EvidenceIntegration.all_alterations_for_detection(
            detection)

        desc_name = ""
        score_dict = {}

        # engineered bring confidence
        if detection.get_engineered() == ENGINEERED:
            score_dict[detection.get_agent()] = detection.get_confidence()
        # natural, turn it blue
        elif detection.get_engineered() == NATURAL:
            score_dict[detection.get_agent()] = detection.get_confidence() * -1
        # indeterminate make it 0
        else:
            score_dict[detection.get_agent()] = -0.25

        if len(all_alterations):
            alteration_group = AlterationGroup(all_alterations[0])
            for alteration in all_alterations[1:]:
                alteration_group.add_alteration(alteration)

            for alteration in all_alterations:
                confidence = -100
                try:
                    confidence = alteration.get_confidence()
                except AttributeError:
                    pass
                score_dict[alteration.get_agent(
                )] = confidence

            try:
                desc_name = desc_name + alteration_group.descriptive_name()
            except AttributeError:
                desc_name = desc_name + str(alteration_group)

        # it has a RegionGroup and AlterationGroup
        if len(all_regions):
            region_group = RegionGroup(all_regions[0])
            for region in all_regions[1:]:
                region_group.add_region(region)

            for region in all_regions:
                confidence = -100
                try:
                    confidence = region.get_confidence()
                except AttributeError:
                    pass
                score_dict[region.get_agent()] = confidence

            # getting scores is always okay for region, but don't add its name to the row
            # label if we already got it from the alteration group.
            if not len(all_alterations):
                try:
                    desc_name = desc_name + region_group.descriptive_name()
                except AttributeError as ae:
                    logging.warn(
                        "Using str(RegionGroup) due to AttributeError: %s", str(ae))
                    desc_name = desc_name + str(region_group)

        return desc_name, score_dict

    @staticmethod
    def detections_to_heatmap_df(all_detections, sample_id=None):
        detections_of_interest = []
        df = pd.DataFrame(
            columns=DetectionsToVisualization.SAMPLE_OVERVIEW_COLUMNS)
        rows = []

        feature_group_detections = []

        for detection in all_detections:
            try:
                if detection.get_agent() == EvidenceIntegration.GUARDIAN_FEATURE_GROUP_AGENT_TAG:
                    feature_group_detections.append(detection)
                    continue

                if sample_id is not None:
                    if detection.get_sample() != sample_id:
                        continue

                if detection.get_singular() != True:
                    continue

                detections_of_interest.append(detection)
            except AttributeError as ae:
                logging.warn("%s, continuing.", str(ae))
                continue

        # list of dicts mapping agent->confidence
        data = []
        # list of row labels - the short name for what was detected
        index = []
        for detection in detections_of_interest:
            desc_name, score_dict = DetectionsToVisualization.name_and_scores_for_composite_detection(
                detection)
            # logging.debug("desc_name: %s", desc_name)
            # logging.debug("score_dict: %s", str(score_dict))
            data.append(score_dict)
            # TODO BBAS: remove this first 20 character limiter once desc_name returned is actually good
            index.append(desc_name)

        df = pd.DataFrame(data, index=index)
        return df
