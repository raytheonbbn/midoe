import copy
import logging
import os
import sys
import json
from schema.dataschema import Detection, Investigation, dumpJSON, jsonobjproc, Feature


class PreprocessInvestigations():

    @staticmethod
    def process_investigations(all_investigations):
        processed_investigations = [
            copy.deepcopy(i) for i in all_investigations]
        for investigation in processed_investigations:
            detections = investigation.get_detections()
            for detection in detections:
                PreprocessInvestigations.cascade_agent_detection(detection)

        PreprocessInvestigations.assign_agent_union_features(
            processed_investigations)

        # one contig per detection at this time... 
        c_split_investigations = PreprocessInvestigations.split_detection_contigs(processed_investigations)


        # one alteration per detection at this time...
        a_split_investigations = PreprocessInvestigations.split_detection_alterations(processed_investigations)

        # NOT DOING READS AT THE MOMENT... TARGETED SEARCH INTRODUCES TOO MANY... COME BACK LATER AND INVESTIGATE
        # THE IMPLICATIONS.... VERY UNEASY

        return c_split_investigations + a_split_investigations


    @staticmethod
    def split_detection_alterations(processed_investigations):
        a_split_investigations = [
            copy.deepcopy(i) for i in processed_investigations]

        for investigation in a_split_investigations:
            detections = investigation.get_detections()
            split_detections = []
            for detection in detections:
                for alteration in detection.get_alterations():
                    if PreprocessInvestigations.alteration_is_empty(alteration):
                        continue
                    else:
                        detection_copy = copy.deepcopy(detection)
                        detection_copy.set_alterations([alteration])
                        split_detections.append(detection_copy)

            investigation.set_detections(split_detections)

        return a_split_investigations


    @staticmethod
    def split_detection_contigs(processed_investigations):
        c_split_investigations = [
            copy.deepcopy(i) for i in processed_investigations]

        for investigation in c_split_investigations:
            detections = investigation.get_detections()
            split_detections = []
            for detection in detections:
                for contig in detection.get_contigs():
                    if contig.has_regions():
                        for region in contig.get_regions():
                            # if PreprocessInvestigations.contig_is_empty(contig):
                                # continue
                            # else:
                            detection_copy = copy.deepcopy(detection)
                            contig_copy = copy.deepcopy(contig)
                            contig_copy.set_regions([region])
                            detection_copy.set_contigs([contig_copy])
                            split_detections.append(detection_copy)
                    else:
                        detection_copy = copy.deepcopy(detection)
                        detection_copy.set_contigs([contig])
                        split_detections.append(detection_copy)
                    
            investigation.set_detections(split_detections)

        return c_split_investigations

    @staticmethod
    def alteration_is_empty(alteration):
        # return true if the alteration has any of the 3 following fields set
        type_exists = True
        try:
            alteration.get_type()
        except AttributeError:
            type_exists = False

        derived_feature = True
        try:
            alteration.get_derivedFeature()
        except AttributeError:
            derived_feature = False

        target_assembly = True
        try:
            alteration.get_targetAssembly()
        except AttributeError:
            target_assembly = False


        return not (type_exists or derived_feature or target_assembly)


    @staticmethod
    def contig_is_empty(contig):
        # return true if the contig has nothing but source set
        reads = True
        try:
            contig.get_reads()
        except AttributeError:
            reads = False

        regions = True
        try:
            contig.get_regions()
        except AttributeError:
            regions = False

        engineered = True
        try:
            contig.get_engineered()
        except AttributeError:
            engineered = False


        return not (reads or regions or engineered)


    @staticmethod
    def cascade_agent_detection(detection):
        if detection.get_agent() is None or detection.get_agent() == "":
            raise AttributeError(
                "Can not cascade the agent property down to proper objects if it is not set in the parent Detection.")

        agent_in_detection = detection.get_agent()

        # go down contigs chain
        all_contigs = []
        try:
            if detection.get_contigs() is not None:
                all_contigs = detection.get_contigs()
        except AttributeError:
            pass
        for contig in all_contigs:
            PreprocessInvestigations.cascade_agent_contig(
                contig, agent_in_detection)

        # go down reads chain
        all_reads = []
        try:
            if detection.get_reads() is not None:
                all_reads = detection.get_reads()
        except AttributeError:
            pass
        for read in all_reads:
            PreprocessInvestigations.cascade_agent_read(
                read, agent_in_detection)

        # go down alterations chain
        alterations = []
        try:
            if detection.get_alterations() is not None:
                alterations = detection.get_alterations()
        except:
            pass
        for alteration in alterations:
            PreprocessInvestigations.cascade_agent_alteration(
                alteration, agent_in_detection)

        # go down features
        features = []
        try:
            if detection.get_features() is not None:
                features = detection.get_features()
        except:
            pass
        for feature in features:
            PreprocessInvestigations.cascade_agent_feature(
                feature, agent_in_detection)

        return detection

    def cascade_agent_read(read, agent_in_detection):
        if read.get_regions() is not None:
            regions = read.get_regions()
            for region in regions:
                PreprocessInvestigations.cascade_agent_region(
                    region, agent_in_detection)

    def cascade_agent_contig(contig, agent_in_detection):
        regions = []
        try:
            if contig.get_regions() is not None:
                regions = contig.get_regions()
        except AttributeError:
            pass
        for region in regions:
            PreprocessInvestigations.cascade_agent_region(
                region, agent_in_detection)

    def cascade_agent_feature(feature, agent_in_detection):
        try:
            if not isinstance(feature, Feature):
                return
            if feature.get_agent() is not None and feature.get_agent() != "":
                return
        except AttributeError:
            feature.set_agent(agent_in_detection)

    def cascade_agent_region(region, agent_in_detection):
        try:
            if region.get_agent() is not None and region.get_agent() != "":
                pass
        except AttributeError:
            region.set_agent(agent_in_detection)

        try:
            if region.get_sourceAlteration() is not None:
                alteration = region.get_sourceAlteration()
                PreprocessInvestigations.cascade_agent_alteration(
                    alteration, agent_in_detection)
        except AttributeError:
            pass

        try:
            if region.get_feature() is not None:
                feature = region.get_feature()
                PreprocessInvestigations.cascade_agent_feature(
                    feature, agent_in_detection)
        except:
            pass

    def cascade_agent_alteration(alteration, agent_in_detection):
        try:
            if alteration.get_agent() is not None and alteration.get_agent() != "":
                return
        except AttributeError:
            alteration.set_agent(agent_in_detection)

        # not currently going into targetRegion

        try:
            feature = alteration.get_derivedFeature()
            PreprocessInvestigations.cascade_agent_feature(
                feature, agent_in_detection)
        except AttributeError:
            pass

        try:
            feature = alteration.get_targetFeature()
            PreprocessInvestigations.cascade_agent_feature(
                feature, agent_in_detection)
        except AttributeError:
            pass

    @staticmethod
    def assign_agent_union_features(investigations):
        features_to_return = []
        for investigation in investigations:

            # get the agent from the first detection
            agent = None
            # warn if the agent is not the same for all detections in this investigation
            for detection in investigation.get_detections():
                if agent is None:
                    agent = detection.get_agent()
                    continue
                if agent != detection.get_agent():
                    logging.warn("WARNING: detections in investigation %s belong to more than one agent, using the first agent found: %s.",
                                 str(investigation), agent)

            if agent is not None:
                for feature in investigation.get_features():
                    if feature.get_agent() is None:
                        feature.set_agent(agent)

                    features_to_return.append(feature)

        return features_to_return
