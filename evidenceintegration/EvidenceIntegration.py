import logging
import copy

import os
import argparse
import sys
import json
from datetime import date
from evidenceintegration.utils import EvidenceIntegrationUtils
from evidenceintegration.PreprocessInvestigations import PreprocessInvestigations
from schema.dataschema import Detection, Investigation, dumpJSON, jsonobjproc


class EvidenceIntegration():

    AGENT_FEATURE_PREFIX_SEPARATOR = "_"
    GUARDIAN_AGENT_TAG = "GUARDIAN"
    GUARDIAN_FEATURE_GROUP_AGENT_TAG = "GUARDIAN_FEATURE_GROUP"
    UNKNOWN_SAMPLE_TAG = "UNKNOWN_SAMPLE"

    @staticmethod
    def form_region_groups(contig_groups, read_groups, feature_groups):
        region_groups = []
        regions_to_group = []

        # we are going to make region_groups based on condition 2, then purge the ones that are len == 1
        temp_region_groups = []
        for c_g in contig_groups:

            # Each region belongs to a contig in the same contig group.
            # Also, each region has featureStart and featureEnd properties and the regions specified by those properties are sufficiently overlapped. Alternatively, if one or both regions are missing featureStart and featureEnd properties, then the regions are sufficiently overlapped based on their start and end properties.
            if len(c_g.get_contigs()) > 1:
                all_regions_in_contig_group = []

                for contig in c_g.get_contigs():
                    all_regions_in_contig_group = all_regions_in_contig_group + contig.get_regions()

                for region in all_regions_in_contig_group:
                    added = False

                    for region_group in temp_region_groups:
                        if region_group.overlaps(region):
                            region_group.add_region(region)
                            added = True
                            break

                    if not added:
                        temp_region_groups.append(RegionGroup(region))

            else:
                # will only ever loop once...
                for contig in c_g.get_contigs():
                    try:
                        regions_to_group = regions_to_group + contig.get_regions()
                    except AttributeError:
                        pass
                    

        for t_r_g in temp_region_groups:
            if len(t_r_g.get_regions()) > 1:
                region_groups.append(t_r_g)
            else:
                # means it matched with no other regions in the contig group
                regions_to_group = regions_to_group + t_r_g.get_regions()

        # we are going to make region_groups based on condition 2, then purge the ones that are len == 1
        # this time, based on read_groups, not contig_groups
        temp_region_groups = []
        for r_g in read_groups:

            # Each region belongs to a contig in the same contig group.
            # Also, each region has featureStart and featureEnd properties and the regions specified by those properties are sufficiently overlapped. Alternatively, if one or both regions are missing featureStart and featureEnd properties, then the regions are sufficiently overlapped based on their start and end properties.
            if len(r_g.get_reads()) > 1:
                all_regions_in_read_group = []
                for read in r_g.get_reads():
                    all_regions_in_read_group = all_regions_in_read_group + read.get_regions()

                for region in all_regions_in_read_group:
                    added = False

                    for region_group in temp_region_groups:
                        if region_group.overlaps(region):
                            region_group.add_region(region)
                            added = True
                            break

                    if not added:
                        # "No sufficient overlap for this region in ContigGroup, creating new temp_region_group."
                        temp_region_groups.append(RegionGroup(region))

            else:
                # will only ever loop once...
                for read in r_g.get_reads():
                    regions_to_group = regions_to_group + read.get_regions()

        for t_r_g in temp_region_groups:
            if len(t_r_g.get_regions()) > 1:
                region_groups.append(t_r_g)
            else:
                # means it matched with no other regions in the contig group
                regions_to_group = regions_to_group + t_r_g.get_regions()

        for region in regions_to_group:
            added = False
            for region_group in region_groups:
                if region_group.matches(region, feature_groups) and not added:
                    region_group.add_region(region)
                    added = True
                    break
            if not added:
                region_groups.append(RegionGroup(region))

        return region_groups

    @staticmethod
    def all_regions_for_detection(detection):
        all_regions = []

        # get all regions safely
        all_contigs = []
        all_reads = []
        try:
            all_contigs += detection.get_contigs()
        except AttributeError as ae:
            logging.warning(ae)
            pass
        try:
            all_reads += detection.get_reads()
        except AttributeError as ae:
            logging.warning(ae)
            pass

        try:
            for r in all_reads:
                all_regions = all_regions + r.get_regions()
        except AttributeError:
            pass
        try:
            for c in all_contigs:
                all_regions = all_regions + c.get_regions()
        except AttributeError:
            pass
        return all_regions

    @staticmethod
    def all_feature_ids_for_detection(detection, investigation = None):
        all_feature_ids = []
        try:
            all_feature_ids = [f.get_id() for f in detection.get_features()]
        except AttributeError as ae:
            pass

        # go down alterations
        try:
            for alteration in detection.get_alterations():
                all_feature_ids = all_feature_ids + \
                    EvidenceIntegration.all_feature_ids_for_alteration(
                        alteration)
        except AttributeError:
            pass

        # go down contigs
        try:
            for contig in detection.get_contigs():
                all_feature_ids = all_feature_ids + \
                    EvidenceIntegration.all_feature_ids_for_contig(contig)
        except AttributeError as ae:
            pass

        # go down reads
        try:
            for read in detection.get_reads():
                all_feature_ids = all_feature_ids + \
                    EvidenceIntegration.all_feature_ids_for_read(read)
        except AttributeError:
            pass

        # look for subfeatures, and get them
        if investigation is not None:
            subfeatures_to_add = []
            for f in all_feature_ids:

                feature_objs = investigation.get_features()
                for feature_obj in feature_objs:
                    if feature_obj.get_id() == f:
                        if feature_obj.get_subFeatures() is not None:
                            subfeatures_to_add = subfeatures_to_add + feature_obj.get_subFeatures()

            all_feature_ids = all_feature_ids + subfeatures_to_add

        return all_feature_ids

    @staticmethod
    def all_feature_ids_for_contig(contig):
        all_feature_ids = []
        try:
            for region in contig.get_regions():
                all_feature_ids = all_feature_ids + \
                    EvidenceIntegration.all_feature_ids_for_region(region)
        except AttributeError:
            pass

        return all_feature_ids

    @staticmethod
    def all_feature_ids_for_read(read):
        all_feature_ids = []
        try:
            for region in read.get_regions():
                all_feature_ids = all_feature_ids + \
                    EvidenceIntegration.all_feature_ids_for_region(region)
        except AttributeError:
            pass

        return all_feature_ids

    @staticmethod
    def all_feature_ids_for_alteration(alteration):
        all_feature_ids = []

        try:
            all_feature_ids = EvidenceIntegration.all_feature_ids_for_region(
                alteration.get_targetRegion())
        except AttributeError:
            pass

        try:
            all_feature_ids = all_feature_ids + \
                [alteration.get_derivedFeature()]
        except AttributeError:
            pass

        try:
            all_feature_ids = all_feature_ids + \
                [alteration.get_targetFeature()]
        except AttributeError:
            pass

        return all_feature_ids

    @staticmethod
    def all_feature_ids_for_region(region):
        all_feature_ids = []
        try:
            all_feature_ids = [region.get_feature()]
        except AttributeError:
            pass

        return all_feature_ids

    @staticmethod
    def all_alterations_for_detection(detection):
        all_alterations = []
        try:
            all_alterations = detection.get_alterations()
        except AttributeError:
            pass
        return all_alterations

    @staticmethod
    def form_honed_detection_groups(detection_group, region_groups, alteration_groups, feature_groups):
        # max(|honed_detection_groups|) == |region_groups| + |alteration_groups|
        #   if they are perfectly orthogonal
        honed_detection_groups = []

        for detection in detection_group.get_detections():
            all_regions = EvidenceIntegration.all_regions_for_detection(
                detection)
            all_alterations = EvidenceIntegration.all_alterations_for_detection(
                detection)
            all_features = EvidenceIntegration.all_feature_ids_for_detection(
                detection)

            # go through all the regions and alterations in this detection, if any HonedDetectionGroup has
            # a region_group with that region in it,
            # or, an alteration_group with that alteration in it,
            # then add yourself to it (if you aren't already in there)
            # if you don't find a HDG with a region group that contains this region, make a new one
            h_d_g_added_to = None
            feature_id_in_common = None

            # search for region in region_group matches
            added = False
            for region in all_regions:
                for h_d_g in honed_detection_groups:
                    if h_d_g.contains_region(region) and not added:
                        if h_d_g.contains(detection):
                            h_d_g_added_to = h_d_g
                            added = True
                        else:
                            logging.debug(
                                "A pre-existing HonedDetectionGroup has a RegionGroup that contains one of my regions, adding this detection to that HonedDetectionGroup")
                            h_d_g.add_detection(detection)
                            h_d_g_added_to = h_d_g
                            added = True

            # search for alteration in alteration_group matches
            for alteration in all_alterations:
                for h_d_g in honed_detection_groups:
                    if h_d_g.contains_alteration(alteration) and not added:
                        if h_d_g.contains(detection):
                            h_d_g_added_to = h_d_g
                            added = True
                        else:
                            logging.debug(
                                "A pre-existing HonedDetectionGroup has an AlterationGroup that contains one of my alterations, adding this detection to that HonedDetectionGroup")
                            h_d_g.add_detection(detection)
                            h_d_g_added_to = h_d_g
                            added = True

            # one of my Features is in a FeatureGroup with a feature from the HonedDetectionGroup
            for feature in all_features:
                for h_d_g in honed_detection_groups:
                    hdg_features = []
                    for d in h_d_g.get_detections():
                        hdg_features = hdg_features + \
                            EvidenceIntegration.all_feature_ids_for_detection(
                                d)
                    for hdg_feature in hdg_features:
                        for f_g in feature_groups:
                            if f_g.contains(feature) and f_g.contains(hdg_feature):
                                if h_d_g.contains(detection):
                                    h_d_g_added_to = h_d_g
                                    added = True
                                elif h_d_g.contains_feature_id(hdg_feature):
                                    logging.debug(
                                        "A pre-existing HonedDetectionGroup has a Detection with a feature in the same FeatureGroup as one of my Features, adding this detection to that HonedDetectionGroup")
                                    h_d_g.add_detection(detection)
                                    h_d_g_added_to = h_d_g
                                    feature_id_in_common = feature
                                    added = True
                                else:
                                    continue

            # if you found an existing group
                    # but you have a RegionGroup that they don't
                    # or an AlterationGroup that they don't
                    # set it
            if added and h_d_g_added_to is not None:
                region_group_to_add = None
                alteration_group_to_add = None
                for region_group in region_groups:
                    for region in all_regions:
                        if region in region_group.get_regions() and region_group.contains_feature_id(feature_id_in_common):
                            region_group_to_add = region_group
                            break

                for alteration_group in alteration_groups:
                    for alteration in all_alterations:
                        if alteration in alteration_group.get_alterations() and alteration_group.contains_feature_id(feature_id_in_common):
                            alteration_group_to_add = alteration_group
                            break

                if h_d_g_added_to.get_region_group() is None and region_group_to_add is not None:
                    logging.debug(
                        "This detection was added to a pre-existing HonedDetectionGroup, "
                        "and it doesn't have a region_group, but I do. "
                        "Setting it to have my region_group.")
                    h_d_g_added_to.set_region_group(region_group_to_add)

                if h_d_g_added_to.get_alteration_group() is None and alteration_group_to_add is not None:
                    logging.debug(
                        "This detection was added to a pre-existing HonedDetectionGroup, "
                        "and it doesn't have a alteration_group, but I do. "
                        "Setting it to have my alteration_group.")
                    h_d_g_added_to.set_alteration_group(
                        alteration_group_to_add)

            # if you do not match any existing group, find your proper region_group and alteration_group, and create
            if not added:
                region_group_to_add = None
                alteration_group_to_add = None
                for region_group in region_groups:
                    for region in all_regions:
                        if region in region_group.get_regions():
                            region_group_to_add = region_group
                            break

                for alteration_group in alteration_groups:
                    for alteration in all_alterations:
                        if alteration in alteration_group.get_alterations():
                            alteration_group_to_add = alteration_group
                            break
                logging.debug(
                    "Creating a new HonedDetectionGroup for this detection because there is no existing HDG with a RegionGroup or AlterationGroup that contains one of my Regions or Alterations.")
                honed_detection_groups.append(HonedDetectionGroup(
                    detection, region_group=region_group_to_add, alteration_group=alteration_group_to_add))
                added = True

        return honed_detection_groups

    @staticmethod
    def form_contig_or_read_groups(entries, contigs=False):
        groups = []

        for x in entries:
            added = False
            for group in groups:
                if group.source_matches(x):
                    group.add(x)
                    added = True
                    break
            if not added:
                if contigs == True:
                    groups.append(ContigGroup(x))
                else:
                    groups.append(ReadGroup(x))

        return groups

    @staticmethod
    def feature_groups_into_detections(prefixed_feature_groups):
        detections = []
        for p_f_g in prefixed_feature_groups:
            detection_dict = {"agent": EvidenceIntegration.GUARDIAN_FEATURE_GROUP_AGENT_TAG,
                              "date": str(date.today()),
                              "singular": True,
                              "features": p_f_g.get_features()}
            detection = Detection(detection_dict)
            detections.append(detection)
        return detections

    @staticmethod
    def form_alteration_groups(detection_group, feature_groups):
        all_alterations = []
        alteration_groups = []
        for d in detection_group.get_detections():
            try:
                all_alterations = all_alterations + d.get_alterations()
            except AttributeError:
                continue

        for alteration in all_alterations:
            added = False
            for alteration_group in alteration_groups:
                if alteration_group.matches(alteration, feature_groups) and not added:
                    alteration_group.add_alteration(alteration)
                    added = True
                    break
            if not added:
                alteration_groups.append(AlterationGroup(alteration))

        return alteration_groups

    @staticmethod
    def form_feature_groups(features):
        feature_groups = []
        features_to_group = features
        keep_looping = 1
        previous_loop_features_remaining = -1

        while keep_looping:
            features_to_group_next = []
            for feature in features_to_group:
                added = False
                for feature_group in feature_groups:
                    if feature_group.matches(feature) and not added:
                        feature_group.add_feature(feature)
                        added = True
                        break
                if not added:
                    features_to_group_next.append(feature)

            if previous_loop_features_remaining == len(features_to_group_next) and len(features_to_group_next) != 0:
                feature_groups.append(FeatureGroup(features_to_group_next[0]))
                if len(features_to_group_next) == 1:
                    return feature_groups
                else:
                    features_to_group_next = features_to_group_next[1:]

            features_to_group = features_to_group_next
            keep_looping = len(features_to_group_next)
            previous_loop_features_remaining = keep_looping

        return feature_groups

    @staticmethod
    def form_detection_groups(detections):
        detection_groups = []

        for detection in detections:
            added = False

            # first we want to look through each detection group, and see if there are any
            # with identical an identical "sample" property
            for detection_group in detection_groups:
                if detection_group.sample_matches(detection):
                    detection_group.add_detection(detection)
                    added = True
                    break

            # second we want to look through each detection group, and see if there are any
            # with identical fastq properties
            if not added:
                for detection_group in detection_groups:
                    if detection_group.fastq_matches(detection):
                        detection_group.add_detection(detection)
                        added = True
                        break

            # third we want to look through each detection group, and see if there are any
            # with identical fasta properties
            if not added:
                for detection_group in detection_groups:
                    if detection_group.fasta_matches(detection):
                        detection_group.add_detection(detection)
                        added = True
                        break

            if not added:
                detection_groups.append(DetectionGroup(detection))

        return detection_groups

    @staticmethod
    def union_features(investigations):
        all_features = []
        for investigation in investigations:
            all_features += investigation.get_features()
        return all_features

    @staticmethod
    def union_detections(investigations):
        all_detections = []
        for investigation in investigations:
            all_detections += investigation.get_detections()
        return all_detections

    @staticmethod
    def integrate(all_investigations, sample_id=None, json_out_dir=None):
        sample_to_guardian_investigation = {}
        # 0 : detections
        # 1 : features
        sample_to_detections_features = {}

        for investigation in all_investigations:
            for detection in investigation.get_detections():

                if (detection.get_agent() == EvidenceIntegration.GUARDIAN_AGENT_TAG or
                        detection.get_agent() == EvidenceIntegration.GUARDIAN_FEATURE_GROUP_AGENT_TAG):
                    continue

                try:
                    sample = detection.get_sample()
                except AttributeError:
                    sample = EvidenceIntegration.UNKNOWN_SAMPLE_TAG

                if sample_id is not None:
                    if sample != sample_id:
                        continue

                if sample in sample_to_detections_features:
                    sample_to_detections_features[sample][0].append(detection)
                    feature_ids_for_sample = EvidenceIntegration.all_feature_ids_for_detection(
                        detection, investigation)
                    all_features = investigation.get_features()
                    features_for_sample = []
                    for feature in all_features:
                        if feature.get_id() in feature_ids_for_sample:
                            features_for_sample.append(feature)

                    x = sample_to_detections_features[sample][1] + \
                        features_for_sample
                    sample_to_detections_features[sample][1] = x
                else:
                    feature_ids_for_sample = EvidenceIntegration.all_feature_ids_for_detection(
                        detection, investigation)
                    all_features = investigation.get_features()
                    features_for_sample = []
                    for feature in all_features:
                        if feature.get_id() in feature_ids_for_sample:
                            features_for_sample.append(feature)

                    sample_to_detections_features[sample] = [
                        [detection], features_for_sample]


        for sample in sample_to_detections_features:
            all_detections = sample_to_detections_features[sample][0]
            all_features = sample_to_detections_features[sample][1]

            logging.debug("Integrating a new sample: %s", str(sample))

            unique_features = list(set(all_features))

            feature_groups = EvidenceIntegration.form_feature_groups(
                unique_features)
            logging.debug("feature_groups: len=%d: %s", len(
                feature_groups), str(feature_groups))

            detection_groups = EvidenceIntegration.form_detection_groups(
                all_detections)
            logging.debug("detection_groups: len=%d: %s", len(
                detection_groups), str(detection_groups))

            guardian_investigation = Investigation()

            for d_g in detection_groups:
                logging.debug("d_g: %s", str(d_g))

                all_contigs = []
                all_reads = []
                for detection in d_g.get_detections():
                    try:
                        all_contigs += detection.get_contigs()
                    except AttributeError:
                        pass
                    try:
                        all_reads += detection.get_reads()
                    except AttributeError:
                        pass

                logging.debug("len(all_reads)=%d", len(all_reads))
                logging.debug("len(all_contigs)=%d", len(all_contigs))
                contig_groups = EvidenceIntegration.form_contig_or_read_groups(
                    all_contigs, contigs=True)
                read_groups = EvidenceIntegration.form_contig_or_read_groups(
                    all_reads, contigs=False)
                logging.debug("len(contig_groups)=%d", len(contig_groups))
                logging.debug("len(read_groups)=%d", len(read_groups))

                # purely debugging, rm before develop merge
                all_regions = []
                for r in all_reads:
                    try:
                        all_regions = all_regions + r.get_regions()
                    except AttributeError:
                        pass
                for c in all_contigs:
                    try:
                        all_regions = all_regions + c.get_regions()
                    except AttributeError:
                        pass
                logging.debug("len(all_regions)=%d", len(all_regions))
                # purely debugging, rm before develop merge

                region_groups = EvidenceIntegration.form_region_groups(
                    contig_groups, read_groups, feature_groups)
                logging.debug("len(region_groups)=%d", len(region_groups))

                alteration_groups = EvidenceIntegration.form_alteration_groups(
                    d_g, feature_groups)
                logging.debug("len(alteration_groups)=%d",
                              len(alteration_groups))

                honed_detection_groups = EvidenceIntegration.form_honed_detection_groups(
                    d_g, region_groups, alteration_groups, feature_groups)
                logging.debug("len(honed_detection_groups)=%d",
                              len(honed_detection_groups))
                for h_d_g in honed_detection_groups:
                    logging.debug("h_d_g: %s", str(h_d_g))

                guardian_composite_detections = []
                for hdg in honed_detection_groups:
                    # [hdg.meld() for hdg in honed_detection_groups]
                    guardian_composite_detections.append(hdg.meld())
                logging.debug("len(guardian_composite_detections)=%d",
                              len(guardian_composite_detections))

                for gcd in guardian_composite_detections:
                    guardian_investigation.add_detection(gcd)

            for f in unique_features:
                guardian_investigation.add_feature(f)

            feature_group_detections = EvidenceIntegration.feature_groups_into_detections(
                feature_groups)

            for fgd in feature_group_detections:
                features = 0
                try:
                    features = len(fgd.get_features())
                except:
                    pass
                logging.debug(
                    "Adding a Detection for a feature group. with %d features", features)
                guardian_investigation.add_detection(fgd)

            logging.debug("guardian_investigation: %s",
                          str(guardian_investigation))

            sample_to_guardian_investigation[sample] = guardian_investigation

        for sample in sample_to_guardian_investigation:
            if json_out_dir is not None:
                json_out = json_out_dir + "/" + str(sample) + ".json"
                logging.debug("Writing JSON Result for sample %s to: %s", str(
                    sample), str(json_out_dir))
                dumpJSON(
                    sample_to_guardian_investigation[sample], open(json_out, "w"))

        return sample_to_guardian_investigation


class ContigGroup(object):
    def __init__(self, contig):
        self.contigs = [contig]

    def __str__(self):
        return "ContigGroup of len=" + str(len(self.contigs))

    def add(self, contig):
        self.contigs.append(contig)

    def get_contigs(self):
        return self.contigs

    def source_matches(self, contig):
        for c in self.contigs:
            try:
                if c.get_source() == contig.get_source():
                    return True
            except AttributeError as ae:
                logging.debug(
                    "One of the contigs does not have source set.  Continuing")
        return False


class ReadGroup(object):
    def __init__(self, read):
        self.reads = [read]

    def __str__(self):
        return "ReadGroup of len=" + str(len(self.reads))

    def add(self, read):
        self.reads.append(read)

    def get_reads(self):
        return self.reads

    def source_matches(self, read):
        for r in self.reads:
            try:
                if r.get_source() == read.get_source():
                    return True
            except AttributeError as ae:
                pass
        return False


class DetectionGroup(object):
    def __init__(self, detection):
        self.detections = [detection]

    def __str__(self):
        return "DetectionGroup of len=" + str(len(self.detections))

    def add_detection(self, detection):
        self.detections.append(detection)

    def get_detections(self):
        return self.detections

    def contains(self, detection):
        return detection in self.detections

    def fasta_matches(self, detection):
        for d in self.detections:
            try:
                if d.get_fasta() == detection.get_fasta():
                    return True
            except AttributeError as ae:
                pass
        return False

    def fastq_matches(self, detection):
        for d in self.detections:
            try:
                if d.get_fastq() == detection.get_fastq():
                    return True
            except AttributeError as ae:
                pass
        return False

    def sample_matches(self, detection):
        for d in self.detections:
            if d.get_sample() == detection.get_sample():
                return True
        return False


class HonedDetectionGroup(DetectionGroup):
    def __init__(self, detection, region_group=None, alteration_group=None):
        super().__init__(detection)
        self.region_group = region_group
        self.alteration_group = alteration_group
        logging.debug("HonedDetectionGroup created with region_group='%s', and alteration_group='%s'", str(
            region_group), str(alteration_group))

    def __str__(self):
        string = "HonedDetectionGroup consists of %d detections, " % (
            len(self.detections))
        if self.region_group is not None:
            string = string + ("a RegionGroup of len=%d, " %
                               len(self.region_group.get_regions()))
        if self.alteration_group is not None:
            string = string + ("a AlterationGroup of len=%d" %
                               len(self.alteration_group.get_alterations()))
        return string

    def meld(self):
        # lets just take the first detection, make a deep copy, and do what we need from there.
        if len(self.detections) < 1:
            raise ValueError(
                "Can not meld HonedDetectionGroup in to a new detection if there are no base detections.")
        detection_result = copy.deepcopy(self.detections[0])

        # Create a new detection with the agent property “GUARDIAN”.
        # Set the singular property of the new detection to true.
        detection_result.set_agent(EvidenceIntegration.GUARDIAN_AGENT_TAG)
        detection_result.set_singular(True)

        # average the confidences
        confidence_sum = 0
        i = 0
        for d in self.detections:
            try:
                confidence_sum = confidence_sum + d.get_confidence()
                i = i + 1
            except AttributeError:
                continue
        try:
            detection_result.set_confidence(confidence_sum / i)
        except ZeroDivisionError:
            # no confidence -> set to 0
            detection_result.set_confidence(0)

        # Then assign that detection the same sample, fasta, and/or fastq properties as the detections in this detection group.
        # If there are non-identical properties among these detections, then leave that property blank in the new detection.
        for d in self.detections:
            try:
                if detection_result.get_sample() != d.get_sample():
                    detection_result.set_sample("")
                    break
            except AttributeError:
                detection_result.set_sample("")
                break
        for d in self.detections:
            try:
                if detection_result.get_fastq() != d.get_fastq():
                    detection_result.set_fastq("")
                    break
            except AttributeError:
                detection_result.set_fastq("")
                break
        for d in self.detections:
            try:
                if detection_result.get_fasta() != d.get_fasta():
                    detection_result.set_fasta("")
                    break
            except AttributeError:
                detection_result.set_fasta("")
                break

        # Add copies of all regions and their mapped read/contigs in the meta-group’s region group to the
        # new detection’s lists of reads and contigs.
        if self.region_group is not None:
            reads_to_add = []
            contigs_to_add = []
            # if a detection in the HonedDetectionGroup has a read, and that read then has a region in my region_group,
            # we are going to add it to the new detection_result

            # TODO: ensure that the entire region_group gets covered in this fashion... we may need to do keep track of what
            #       read a region came from if this occurs...
            for detection in self.detections:
                detection_reads = []
                try:
                    detection_reads = detection.get_reads()
                except AttributeError:
                    pass

                for read in detection_reads:
                    for region in read.get_regions():
                        if self.contains_region(region) and read not in reads_to_add:
                            reads_to_add.append(read)

                # same thing for contigs
                detection_contigs = []
                try:
                    detection_contigs = detection.get_contigs()
                except AttributeError:
                    pass
                for contig in detection_contigs:
                    regions = []
                    try:
                        regions = contig.get_regions()
                    except AttributeError:
                        pass
                    for region in regions:
                        if self.contains_region(region) and contig not in contigs_to_add:
                            contigs_to_add.append(contig)

            detection_result.set_reads(reads_to_add)
            detection_result.set_contigs(contigs_to_add)

        # Add copies of all alterations in the meta-group’s alteration group to the new detection’s list of alterations.
        if self.alteration_group is not None:
            detection_result.set_alterations(
                self.alteration_group.get_alterations())

        return detection_result

    def contains_region(self, region):
        if self.region_group is None:
            return False
        if region in self.region_group.get_regions():
            return True
        else:
            return False

    def contains_alteration(self, alteration):
        if self.alteration_group is None:
            return False
        if alteration in self.alteration_group.get_alterations():
            return True
        else:
            return False

    def contains_feature_id(self, feature_id):
        if self.alteration_group is not None:
            if self.alteration_group.contains_feature_id(feature_id):
                return True

        if self.region_group is not None:
            if self.region_group.contains_feature_id(feature_id):
                return True

        return False

    def get_alteration_group(self):
        return self.alteration_group

    def get_region_group(self):
        return self.region_group

    def set_alteration_group(self, alteration_group):
        if self.alteration_group is not None:
            raise Exception(
                "Unable to re-set AlterationGroup in a HonedDetectionGroup.")
        self.alteration_group = alteration_group

    def set_region_group(self, region_group):
        if self.region_group is not None:
            raise Exception(
                "Unable to re-set RegionGroup in a HonedDetectionGroup.")
        self.region_group = region_group


class AlterationGroup(object):
    def __init__(self, alteration):
        self.alterations = [alteration]

    def __str__(self):
        result = "AlterationGroup, len=" + str(len(self.alterations))
        result = result + ", alterations:["

        # all but the last alteration
        for alteration in self.alterations[:-1]:
            result = result + str(alteration) + "; "
        # the last alteration
        for alteration in self.alterations[-1:]:
            result = result + str(alteration) + "]"

        return result

    def descriptive_name(self):
        result = ""

        # for now, we are going to get the first alteration, and create the name from it...
        # if the alterations converge into an alteration group they should match descriptive
        # names.
        alteration = self.alterations[0]

        result = result + alteration.get_type()

        if alteration.get_type().lower() == "insertion":
            result = result + " of "
            result = result + alteration.get_derivedFeature()

        if alteration.get_type().lower() == "deletion":
            result = result + " in "
            result = result + alteration.get_targetFeature()

        return result

    def descriptive_name_for_csv(self):
        result = ""
        alteration = self.alterations[0]

        if alteration.get_type().lower() == "insertion":
            shortest_feature = min([alteration.get_derivedFeature() for alteration in self.alterations], key=len)

            for alteration in self.alterations:
                if alteration.get_derivedFeature() == shortest_feature:
                    if alteration.get_agent() is not None:
                        result = "[" + alteration.get_agent() + "]"
                    result += alteration.get_type() + " of " + alteration.get_derivedFeature()




        if alteration.get_type().lower() == "deletion":
            shortest_feature = min([alteration.get_targetFeature() for alteration in self.alterations], key=len)

            for alteration in self.alterations:
                if alteration.get_targetFeature() == shortest_feature:
                    if alteration.get_agent() is not None:
                        result = "[" + alteration.get_agent() + "]"
                    result += alteration.get_type() + " in " + alteration.get_targetFeature()




        return result


    def add_alteration(self, alteration):
        self.alterations.append(alteration)

    def get_alterations(self):
        return self.alterations

    def contains_feature_id(self, feature_id):
        for alteration in self.alterations:
            try:
                if alteration.get_targetRegion().get_feature() == feature_id:
                    return True
            except:
                pass

            try:
                if alteration.get_derivedFeature() == feature_id:
                    return True
            except:
                pass

            try:
                if alteration.get_targetFeature() == feature_id:
                    return True
            except:
                pass

        return False

    def matches(self, alteration, feature_groups):
        for a in self.alterations:
            if a.rough_equality(alteration, feature_groups):
                return True
        return False


class RegionGroup(object):
    def __init__(self, region):
        self.regions = [region]

    def __str__(self):
        result = "RegionGroup, len=" + str(len(self.regions))
        result = result + ", regions:["

        # all but the last region
        for region in self.regions[:-1]:
            result = result + str(region) + "; "
        # the last region
        for region in self.regions[-1:]:
            result = result + str(region) + "]"

        return result

    def descriptive_name(self):
        result = "engineered sequence of " + self.regions[0].get_feature()
        return result

    def descriptive_name_for_csv(self):
        result = ""
        shortest_feature = min([region.get_feature() for region in self.regions], key=len)

        for region in self.regions:
            if region.get_feature() == shortest_feature:
                if region.get_agent() is not None:
                    result = "[" + region.get_agent() + "]"
                result += region.get_feature()

        return result

    def add_region(self, region):
        self.regions.append(region)

    def get_regions(self):
        return self.regions

    def overlaps(self, region):
        for r in self.regions:
            if r.overlaps(region):
                return True
        return False

    def contains_feature_id(self, feature_id):
        for region in self.regions:
            try:
                if region.get_feature() == feature_id:
                    return True
            except AttributeError:
                return False
        return False

    def matches(self, region, feature_groups):
        for r in self.regions:
            if r.rough_equality(region, feature_groups):
                return True
        return False


class FeatureGroup(object):
    def __init__(self, feature):
        self.features = [feature]

    def __str__(self):
        return "FeatureGroup of len=" + str(len(self.features))

    def contains(self, feature_id):
        for feature in self.features:
            if feature_id == feature.get_id():
                return True
        return False

    def get_ids(self):
        return [f.get_id() for f in self.features]

    def add_feature(self, feature):
        self.features.append(feature)

    def get_features(self):
        return self.features

    def matches(self, feature):
        for f in self.features:
            if f.rough_equality(feature):
                return True
        return False


def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser()
    parser.add_argument('-j', '--json_files', nargs='*', default=[])
    parser.add_argument('-d', '--json_dirs', nargs='*', default=[])
    parser.add_argument('-s', '--sample_id', nargs=1, default=None)
    parser.add_argument('-l', '--integration_log', nargs='?', default='')

    args = parser.parse_args(args)
    json_files = args.json_files
    json_dirs = args.json_dirs
    sample_id = args.sample_id[0]
    json_out_dir = None

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
                    json_dir, sample_id)
            json_out_dir = json_dir + "/out/"

    for j in json_files:
        if "guardian_evidence.json" in j:
            json_files.remove(j)

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

    processed_investigations = PreprocessInvestigations.process_investigations(
        all_investigations)

    logging.debug("len processed_investigations: %s", str(len(processed_investigations)))

    sample_to_guardian_investigation = EvidenceIntegration.integrate(
        processed_investigations, json_out_dir=json_out_dir)

    for sample in sample_to_guardian_investigation:
        logging.debug("\tsample: %s", str(sample))
        logging.debug("\tguardian_investigation: %s", str(
            sample_to_guardian_investigation[sample]))

    print('Finished')


if __name__ == '__main__':
    main()
