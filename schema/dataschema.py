import logging
import os
import sys
import numpy as np
import json
from Bio.Seq import Seq
from Bio import Align

ENGINEERING_VALS = ["http://guardian.bbn.technology#natural",
                    "http://guardian.bbn.technology#indeterminate_origin",
                    "http://guardian.bbn.technology#engineered"]

ORIENTATION_VALS = ["http://sbols.org/v2#inline",
                    "http://sbols.org/v2#reverseComplement"]


class Investigation(object):
    def __init__(self):
        self.detections = []  # list of Detection objects
        self.features = []  # list of Feature objects

    def __str__(self):
        return "Investigation object with %d detections and %d features." % (len(self.detections), len(self.features))

    def add_detection(self, det):
        self.detections.append(det)

    def add_feature(self, feat):
        self.features.append(feat)

    def read_from_json(self, json_dict):
        if "detections" in json_dict:
            for detection_dict in json_dict["detections"]:
                self.add_detection(Detection(detection_dict))

        if "features" in json_dict:
            for feat_dict in json_dict["features"]:
                self.add_feature(Feature(feat_dict))

    def set_detections(self, detections):
        self.detections = detections

    def get_detections(self):
        return self.detections

    def get_features(self):
        return self.features


class Evidence(object):
    def __init__(self, evidence_dict):
        for key in evidence_dict:
            if key == "engineered":
                self.set_engineered(evidence_dict[key])
            if key == "confidence":
                self.set_confidence(evidence_dict[key])
            if key == "taxa":
                self.set_taxa(evidence_dict[key])
            if key == "hostTaxa":
                self.set_hostTaxa(evidence_dict[key])
            if key == "description":
                self.set_description(evidence_dict[key])
            if key == "scores":
                self.set_scores(evidence_dict[key])

    def __str__(self):
        result = "Evidence: "
        try:
            result = result + ("engineered=%s" % str(self.engineered))
        except AttributeError:
            pass

        try:
            result = result + (",confidence=%s" % str(self.confidence))
        except AttributeError:
            pass

        try:
            result = result + (",taxa=%s" % str(self.taxa))
        except AttributeError:
            pass

        try:
            result = result + (",hostTaxa=%s" % str(self.hostTaxa))
        except AttributeError:
            pass

        try:
            result = result + (",description=%s" % str(self.description))
        except AttributeError:
            pass

        try:
            if len(self.scores > 0):
                # all but the last score
                result = result + ", scores:["
                for score in self.scores[:-1]:
                    result = result + str(score) + "; "
                # the last score
                for score in self.scores[-1:]:
                    result = result + str(score) + "]"
        except AttributeError:
            pass

        return result

    def set_engineered(self, val):
        # See ENGINEERING_VALS for permitted values
        if val in ENGINEERING_VALS:
            self.engineered = val
        else:
            sys.exit('Bad engineered value')

    def set_confidence(self, val):
        #real [0, 1]
        if val >= 0 and val <= 1:
            self.confidence = val
        else:
            sys.exit('Bad confidence value')

    def set_taxa(self, taxal):
        # opt list of URI refs for taxn classn
        self.taxa = taxal

    def set_hostTaxa(self, taxal):
        # opt list of URI refs for taxn classn
        self.hostTaxa = taxal

    def set_role(self, role):
        # opt string class for detected component (promoter, CDS, deletion, ...)
        self.role = role

    def set_description(self, notes):
        # opt string of notes for evidence
        self.description = notes

    def set_scores(self, scores):
        # scores [0..*]
        self.scores = []
        for score_dict in scores:
            self.scores.append(Score(score_dict))

    def get_engineered(self):
        return self.engineered

    def get_hostTaxa(self):
        return self.hostTaxa

    def get_confidence(self):
        return self.confidence


class Score(object):
    def __init__(self, score_dict):
        # must have value, all others optional [0..1]
        self.set_value(score_dict["value"])
        for key in score_dict:
            if key == "type":
                self.set_type(score_dict[key])
            if key == "minValue":
                self.set_minValue(score_dict[key])
            if key == "maxValue":
                self.set_maxValue(score_dict[key])
            if key == "lowerThreshold":
                self.set_lowerThreshold(score_dict[key])
            if key == "upperThreshold":
                self.set_upperThreshold(score_dict[key])

    def __str__(self):
        result = "Score:"
        try:
            result = result + ("type=%s" % str(self.type))
        except AttributeError:
            pass

        try:
            result = result + (",minValue=%s" % str(self.minValue))
        except AttributeError:
            pass

        try:
            result = result + (",maxValue=%s" % str(self.maxValue))
        except AttributeError:
            pass

        try:
            result = result + (",lowerThreshold=%s" % str(self.lowerThreshold))
        except AttributeError:
            pass

        try:
            result = result + (",upperThreshold=%s" % str(self.upperThreshold))
        except AttributeError:
            pass

        return result

    def set_type(self, val):
        # URI [0..1]
        self.type = value

    def set_minValue(self, val):
        # Number [0..1]
        self.minValue = val

    def set_maxValue(self, val):
        # Number [0..1]
        self.maxValue = val

    def set_lowerThreshold(self, val):
        # Number [0..1]
        self.lowerThreshold = val

    def set_upperThreshold(self, val):
        # Number [0..1]
        self.upperThreshold = val


class Detection(Evidence):
    def __init__(self, detection_dict):
        super().__init__(detection_dict)
        self.agent = detection_dict["agent"]  # URI for "detection" agent
        self.batch = "NOBATCH"
        self.reads = []
        self.alterations = []
        self.contigs = []

        for key in detection_dict:
            if key == "fasta":
                self.set_fasta(detection_dict[key])
            if key == "gff":
                self.set_gff(detection_dict[key])
            if key == "fastq":
                self.set_fastq(detection_dict[key])
            if key == "sample":
                self.set_sample(detection_dict[key])
            if key == "runtime":
                self.set_runtime(detection_dict[key])
            if key == "date":
                self.set_date(detection_dict[key])
            if key == "singular":
                self.set_singular(detection_dict[key])
            if key == "reads":
                self.set_reads(detection_dict[key])
            if key == "contigs":
                self.set_contigs(detection_dict[key])
            if key == "alterations":
                self.set_alterations(detection_dict[key])
            if key == "features":
                self.set_features(detection_dict[key])

    def set_batch(self, val):
        self.batch = val

    def set_sample(self, val):
        # optional URI for sample
        self.sample = val

    def set_fasta(self, val):
        # opt URI file ref
        self.fasta = val

    def set_fastq(self, val):
        # opt list of URI file ref
        self.fastq = val

    def set_gff(self, val):
        # opt URI file ref
        self.gff = val

    def set_runtime(self, val):
        # minutes taken by agent
        if val >= 0:
            self.runtime = val
        else:
            sys.exit('Bad runtime value')

    def inc_runtime(self, val):
        if val >= 0:
            self.runtime += val
        else:
            sys.exit('Bad runtime value')

    def set_date(self, val):
        # date string
        self.date = val

    def set_reads(self, reads):
        read_objs = []
        for read in reads:
            if isinstance(read, ReadReference):
                read_objs.append(read)
            else:
                read_objs.append(ReadReference(read))

        # list of ReadReference objects
        self.reads = read_objs

    def set_singular(self, val):
        self.singular = val

    def set_agent(self, val):
        self.agent = val

    def set_contigs(self, contigs):
        contig_objs = []

        for contig in contigs:
            if isinstance(contig, ContigReference):
                contig_objs.append(contig)
            else:
                contig_objs.append(ContigReference(contig))

        # list of ContigReference objects
        self.contigs = contig_objs

    def set_alterations(self, alterations):
        self.alterations = []
        # list of Alteration objects
        for alteration in alterations:
            if isinstance(alteration, Alteration):
                self.alterations.append(alteration)
            else:
                self.alterations.append(Alteration(alteration))

    def set_features(self, features):
        self.features = []
        # list of feature objects
        for feature in features:
            if isinstance(feature, Feature):
                self.features.append(feature)
            else:
                self.features.append(Feature(feature))

    def get_batch(self):
        return self.batch

    def get_features(self):
        return self.features

    def get_alterations(self):
        # list of Alteration objects
        return self.alterations

    def get_reads(self):
        # list of ReadReference objects
        return self.reads

    def get_contigs(self):
        # list of ContigReference objects
        return self.contigs

    def get_agent(self):
        return self.agent

    def get_sample(self):
        return self.sample

    def get_fastq(self):
        return self.fastq

    def get_fasta(self):
        return self.fasta

    def get_date(self):
        return self.date

    def get_singular(self):
        return self.singular


class ReadReference(Evidence):
    def __init__(self, read_dict):
        super().__init__(read_dict)
        self.source = read_dict["source"]  # URI string for SeqID
        for key in read_dict:
            if key == "regions":
                self.set_regions(read_dict[key])

    def set_regions(self, region_list):
        self.regions = []  # list of Region objects
        try:
            for region_dict in region_list:
                self.regions.append(Region(region_dict))
        except:
            for region_obj in region_list:
                if isinstance(region_obj, Region):
                    self.regions.append(region_obj)

    def get_regions(self):
        return self.regions

    def get_source(self):
        return self.source


class ContigReference(Evidence):
    def __init__(self, contig_dict):
        super().__init__(contig_dict)
        self.source = contig_dict["source"]  # URI string for SeqID
        for key in contig_dict:
            if key == "reads":
                self.set_reads(contig_dict[key])
            if key == "regions":
                self.set_regions(contig_dict[key])

    def set_regions(self, region_list):
        self.regions = []  # list of Region objects
        for region in region_list:
            if isinstance(region, Region):
                self.regions.append(region)
            else:
                self.regions.append(Region(region))

    def set_reads(self, reads):
        # opt list of URIs used in assembly
        self.reads = reads

    def get_source(self):
        return self.source

    def get_regions(self):
        return self.regions

    def has_regions(self):
        return hasattr(self, 'regions') and len(self.regions) > 0


class Region(Evidence):
    def __init__(self, region_dict):
        super().__init__(region_dict)
        self.start = region_dict["start"]  # integer > 0
        self.end = region_dict["end"]  # integer >= start
        self.feature = None
        for key in region_dict:
            if key == "source":
                self.set_source(region_dict[key])
            if key == "orientation":
                self.set_orientation(region_dict[key])
            if key == "role":
                self.set_role(region_dict[key])
            if key == "featureStart":
                self.set_featureStart(region_dict[key])
            if key == "featureEnd":
                self.set_featureEnd(region_dict[key])
            if key == "featureOrientation":
                self.set_featureOrientation(region_dict[key])
            if key == "agent":
                self.set_agent(region_dict[key])
            if key == "feature":
                # string, not object for now
                self.set_feature(region_dict[key])
            if key == "sourceAlteration":
                self.set_sourceAlteration(region_dict[key])

    # exact equality across all fields
    def __eq__(self, x):
        try:
            return (isinstance(x, self.__class__) and
                    self.get_feature() == x.get_feature() and
                    self.get_orientation() == x.get_orientation() and
                    self.get_source() == x.get_source() and
                    self.get_role() == x.get_role() and
                    self.get_featureStart() == x.get_featureStart() and
                    self.get_featureEnd() == x.get_featureEnd() and
                    self.get_featureOrientation() == x.get_featureOrientation() and
                    self.get_sourceAlteration() == x.get_sourceAlteration() and
                    self.get_agent() == x.get_agent())
        except AttributeError:
            return False

    def __str__(self):
        result = "Region: (" + super().__str__() + "),"

        try:
            result = result + ("start=%s" % str(self.start))
        except AttributeError:
            pass

        try:
            result = result + (",end=%s" % str(self.end))
        except AttributeError:
            pass

        try:
            result = result + (",source=%s" % str(self.source))
        except AttributeError:
            pass

        try:
            result = result + (",orientation=%s" % str(self.orientation))
        except AttributeError:
            pass

        try:
            result = result + (",role=%s" % str(self.role))
        except AttributeError:
            pass

        try:
            result = result + (",featureStart=%s" % str(self.featureStart))
        except AttributeError:
            pass

        try:
            result = result + (",featureEnd=%s" % str(self.featureEnd))
        except AttributeError:
            pass

        try:
            result = result + (",featureOrientation=%s" %
                               str(self.featureOrientation))
        except AttributeError:
            pass

        try:
            result = result + (",agent=%s" % str(self.agent))
        except AttributeError:
            pass

        try:
            result = result + (",sourceAlteration=%s" %
                               str(self.sourceAlteration))
        except AttributeError:
            pass

        try:
            result = result + (",feature=%s" % str(self.feature))
        except AttributeError:
            pass

        return result

    # Two regions are equivalent if one of the following conditions is satisfied:
        # Each region identifies the same feature.
        # The regions identify features that belong to the same feature group.
        # Each region belongs to a read in the same read group or a contig in the same contig group. Also, each region has featureStart and featureEnd properties and the regions specified by those properties are sufficiently overlapped. Alternatively, if one or both regions are missing featureStart and featureEnd properties, then the regions are sufficiently overlapped based on their start and end properties.
    def rough_equality(self, other_region, feature_groups=None):
        # Each region identifies the same feature.
        condition_0 = False
        try:
            condition_0 = (self.get_feature() == other_region.get_feature())
        except AttributeError:
            # print("Either self or other_region does not have a feature set, condition_0 = False; continuing.")
            pass

        # The regions identify features that belong to the same feature group.
        condition_1 = False
        if feature_groups is not None:
            for f_g in feature_groups:
                try:
                    if f_g.contains(self.get_feature()) and f_g.contains(other_region.get_feature()):
                        condition_1 = True
                except AttributeError as ae:
                    # print("Either self or other_region does not have a feature set; continuing.")
                    continue

        # Each region has a sourceAlteration that identifies the same targetFeature
        #   or a targetFeature that belongs to the same feature group.
        condition_3 = False
        if feature_groups is not None:
            for f_g in feature_groups:
                try:
                    if f_g.contains(self.get_sourceAlteration().get_targetFeature()) and f_g.contains(other_region.get_sourceAlteration().get_targetFeature()):
                        condition_3 = True
                except AttributeError as ae:
                    continue

        return (isinstance(other_region, self.__class__) and
                (condition_0 or condition_1 or condition_3)
                )

    def overlaps(self, other_region):
        try:
            # self covers other_region
            if self.get_featureStart() <= other_region.get_featureStart() and self.get_featureEnd() >= other_region.get_featureEnd():
                return True

            # other_region covers self
            if other_region.get_featureStart() <= self.get_featureStart() and other_region.get_featureEnd() >= self.get_featureEnd():
                return True

            return False
        except AttributeError:
            return False

    def set_orientation(self, val):
        #inline or reverseComplement
        if val in ORIENTATION_VALS:
            self.orientation = val
        else:
            print("orientation val not in ORIENTATION_VALS, not assigning.")
            return

    def set_source(self, val):
        # URI for GFF annotation
        self.source = val

    def set_role(self, val):
        # string (term from sequence_feature branch)
        self.role = val

    def set_feature(self, val):
        # string ID ref to Feature in parent Investigation
        self.feature = val

    def set_featureStart(self, val):
        # int start for associated feature
        self.featureStart = val

    def set_featureEnd(self, val):
        # int end for associated feature
        self.featureEnd = val

    def set_featureOrientation(self, val):
        #inline or reverseComplement
        if val in ORIENTATION_VALS:
            self.featureOrientation = val
        else:
            print("orientation val not in ORIENTATION_VALS, not assigning.")
            return

    def set_sourceAlteration(self, alt_dict):
        # Alteration object
        self.sourceAlteration = Alteration(alt_dict)

    def set_agent(self, val):
        # URI for agent that detected this
        self.agent = val

    def get_start(self):
        return self.start

    def get_end(self):
        return self.end

    def get_feature(self):
        return self.feature

    def get_orientation(self):
        return self.orientation

    def get_source(self):
        return self.source

    def get_role(self):
        return self.role

    def get_featureStart(self):
        return self.featureStart

    def get_featureEnd(self):
        return self.featureEnd

    def get_featureOrientation(self):
        return self.featureOrientation

    def get_sourceAlteration(self):
        return self.sourceAlteration

    def get_agent(self):
        return self.agent


class Feature(object):
    def __init__(self, feat_dict):
        self.id = feat_dict["id"]  # feature ID string
        self.role = None
        self.source = None
        self.sourceStart = None
        self.sourceEnd = None
        self.sourceOrientation = None
        self.sequence = None
        self.name = None
        self.subFeatures = None
        self.agent = None

        for key in feat_dict:
            if key == "name":
                self.set_name(feat_dict[key])  # optional, string, name
            if key == "source":
                # opt URI ref to external record e.g. NCBI SeqID
                self.set_source(feat_dict[key])
            if key == "sourceStart":
                # int start in context of external record
                self.set_sourceStart(feat_dict[key])
            if key == "sourceEnd":
                # int end in context of external record
                self.set_sourceEnd(feat_dict[key])
            if key == "sourceOrientation":
                # inline or reverseComplement
                self.set_sourceOrientation(feat_dict[key])
            if key == "sequence":
                self.set_sequence(feat_dict[key])  # string of IUPAC DNA codes
            if key == "role":
                # string (term from sequence_feature branch)
                self.set_role(feat_dict[key])
            if key == "agent":
                # URI for agent that detected this
                self.set_agent(feat_dict[key])
            if key == "subFeatures":
                self.set_subFeatures(feat_dict[key])

    def __str__(self):
        return "Feature:id=%s;name=%s;source=%s;role=%s;agent=%s" % (
            str(self.id), str(self.name), str(self.source), str(self.role), str(self.agent))

    def __hash__(self):
        return hash(self.id)

    # exact equality across all fields
    def __eq__(self, x):
        return (isinstance(x, self.__class__) and
                self.get_role() == x.get_role() and
                self.get_source() == x.get_source() and
                self.get_sourceStart() == x.get_sourceStart() and
                self.get_sourceEnd() == x.get_sourceEnd() and
                self.get_sourceOrientation() == x.get_sourceOrientation() and
                self.get_sequence() == x.get_sequence() and
                self.get_name() == x.get_name() and
                self.get_subFeatures() == x.get_subFeatures() and
                self.get_agent() == x.get_agent())

    # rough equality, defined as:
        # Each feature has a sequence and those sequences are sufficiently similar.
        # Each feature has a source property and those source properties are identical, and each feature has sourceStart and sourceEnd properties and the regions specified by these properties are sufficiently overlapped.
    def rough_equality(self, other_feature):

        # seq_match = self.sequences_sufficiently_similar(other_feature)
        # if seq_match:
        #     logging.debug("\n\n\n\n")
        # logging.debug("feature '%s' vs feature '%s', sequences_sufficiently_similar=%s", self.get_id(), other_feature.get_id(), seq_match)

        return (isinstance(other_feature, self.__class__) and
                (
                # Each feature has a sequence and those sequences are sufficiently similar.
                (self.sequences_sufficiently_similar(other_feature)) or

                (
                    # Each feature has a source property and those source properties are identical.
                    self.get_source() == other_feature.get_source() and

                    # and each feature has sourceStart and sourceEnd properties and the regions specified by these properties are sufficiently overlapped.
                    self.sources_sufficiently_overlap(other_feature))
                )
                )

    def sequences_sufficiently_similar(self, other_feature):

        if (self.get_sequence() is None or other_feature.get_sequence() is None):
            return False

        seq_0 = self.get_sequence().upper()
        seq_1 = other_feature.get_sequence().upper()

        if len(seq_0) < 20 or len(seq_1) < 20:
            return False
        elif seq_0 == seq_1: # exact sequence match
            return True
        elif seq_0 in seq_1 or seq_1 in seq_0: # either sequence contains the other
            return True
        elif (len(seq_0) < 1000 or len(seq_1) < 1000): # perform fuzzy sequence match
            score_threshold = 1 - min(len(seq_0), len(seq_1))*0.0005
        else:
            score_threshold = 0.5

        rc_seq_0 = str(Seq(seq_0).reverse_complement())

        aligner = Align.PairwiseAligner()
        aligner.match_score = 1
        aligner.mismatch_score = -2
        aligner.internal_gap_score = -2.5

        score = aligner.score(seq_0, seq_1)
        score_scaled_for_len = score / min(len(seq_0), len(seq_1))
        rc_score = aligner.score(rc_seq_0, seq_1)
        rc_score_scaled_for_len = rc_score / min(len(seq_0), len(seq_1))
        # print("score: ", score)
        # print("\t\t\t\tscore_scaled_for_len: ", score_scaled_for_len)
        # print("rc_score: ", rc_score)
        # print("\t\t\t\trc_score_scaled_for_len: ", rc_score_scaled_for_len)

        return score_scaled_for_len > score_threshold or rc_score_scaled_for_len > score_threshold

    def sources_sufficiently_overlap(self, other_feature):
        try:
            # self covers other_feature
            if self.get_sourceStart() <= other_feature.get_sourceStart() and self.get_sourceEnd() >= other_feature.get_sourceEnd():
                return True

            # other_feature covers self
            if other_feature.get_sourceStart() <= self.get_sourceStart() and other_feature.get_sourceEnd() >= self.get_sourceEnd():
                return True

            return False
        except TypeError as te:
            return False

    def set_id(self, val):
        self.id = val

    def set_role(self, val):
        # string (term from sequence_feature branch)
        self.role = val  # string (term from sequence_feature branch)

    def set_source(self, val):
        # opt URI ref to external record e.g. NCBI SeqID
        self.source = val

    def set_sourceStart(self, val):
        # int start in context of external record
        self.sourceStart = val

    def set_sourceEnd(self, val):
        # int end in context of external record
        self.sourceEnd = val

    def set_sourceOrientation(self, val):
        #inline or reverseComplement
        if val in ORIENTATION_VALS:
            self.sourceOrientation = val
        else:
            print("orientation val not in ORIENTATION_VALS, not assigning.")
            return

    def set_sequence(self, val):
        # string of IUPAC DNA codes
        self.sequence = val

    def set_name(self, name):
        # optional string name
        self.name = name

    def set_subFeatures(self, subfeats):
        # optional list of subfeatures
        self.subFeatures = subfeats

    def set_agent(self, val):
        # URI for agent that detected this
        self.agent = val

    def get_id(self):
        return self.id

    def get_role(self):
        return self.role

    def get_source(self):
        return self.source

    def get_sourceStart(self):
        return self.sourceStart

    def get_sourceEnd(self):
        return self.sourceEnd

    def get_sourceOrientation(self):
        return self.sourceOrientation

    def get_sequence(self):
        return self.sequence

    def get_name(self):
        return self.name

    def get_subFeatures(self):
        return self.subFeatures

    def get_agent(self):
        return self.agent


class Alteration(Evidence):
    def __init__(self, alt_dict):
        super().__init__(alt_dict)
        self.type = alt_dict["type"]
        self.derivedFeature = None
        self.targetFeature = None
        self.readCount = None
        for key in alt_dict:
            if key == "targetSeq":
                self.set_targetSeq(alt_dict[key])
            if key == "targetAssembly":
                self.set_targetAssembly(alt_dict[key])
            if key == "agent":
                self.set_agent(alt_dict[key])
            if key == "derivedFeature":
                self.set_derivedFeature(alt_dict[key])
            if key == "targetFeature":
                self.set_targetFeature(alt_dict[key])
            if key == "readCount":
                self.set_readCount(alt_dict[key])

    # exact equality across all fields
    # TODO, parent class __eq__, call it here and every class that inherits Evidence
    def __eq__(self, x):
        try:
            return (isinstance(x, self.__class__) and
                    self.get_type() == x.get_type() and
                    self.get_derivedFeature() == x.get_derivedFeature() and
                    self.get_targetSeq() == x.get_targetSeq() and
                    self.get_targetAssembly() == x.get_targetAssembly() and
                    self.get_targetFeature() == x.get_targetFeature() and
                    self.get_targetRegion() == x.get_targetRegion() and
                    self.get_agent() == x.get_agent())
        except AttributeError:
            return False

    # rough equality, defined as:
        # must have identical types and one of the following conditions is satisfied:
        # 0 The alterations identify the same derivedFeature
        # 1 they identify derivedFeatures that belong to the same feature group.
        # 2 they identify the same targetFeature
        # 3 they identify targetFeatures that belong to the same feature group.
    def rough_equality(self, other_alteration, feature_groups):

        # The alterations identify the same derivedFeature
        condition_0 = False
        try:
            if self.get_derivedFeature() == other_alteration.get_derivedFeature():
                condition_0 = True
        except AttributeError:
            condition_0 = False

        # 1 they identify derivedFeatures that belong to the same feature group.
        condition_1 = False
        for f_g in feature_groups:
            try:
                if f_g.contains(self.get_derivedFeature()) and f_g.contains(other_alteration.get_derivedFeature()):
                    condition_1 = True
            except AttributeError:
                condition_1 = False

        # they identify the same targetFeature
        condition_2 = False
        try:
            if (self.get_targetFeature() is not None) and self.get_targetFeature() == other_alteration.get_targetFeature():
                condition_2 = True
        except AttributeError:
            condition_2 = False

        # 3 they identify targetFeatures that belong to the same feature group.
        condition_3 = False
        for f_g in feature_groups:
            try:
                if f_g.contains(self.get_targetFeature()) and f_g.contains(other_alteration.get_targetFeature()):
                    condition_3 = True
            except AttributeError:
                condition_3 = False

        return (isinstance(other_alteration, self.__class__) and
                (self.get_type() == other_alteration.get_type()) and
                (
                    # The alterations identify the same derivedFeature
                    condition_0 or

                    # they identify derivedFeatures that belong to the same feature group.
                    condition_1 or

                    # they identify the same targetFeature
                    condition_2 or

                    # they identify targetFeatures that belong to the same feature group.
                    condition_3

        )
        )

    def __str__(self):
        result = "Alteration: (" + super().__str__() + "),"

        try:
            result = result + ("type=%s" % str(self.type))
        except AttributeError:
            pass

        try:
            result = result + (",targetSeq=%s" % str(self.targetSeq))
        except AttributeError:
            pass

        try:
            result = result + (",targetAssembly=%s" % str(self.targetAssembly))
        except AttributeError:
            pass

        try:
            result = result + (",agent=%s" % str(self.agent))
        except AttributeError:
            pass

        try:
            result = result + (",derivedFeature=%s" % str(self.derivedFeature))
        except AttributeError:
            pass

        try:
            result = result + (",targetFeature=%s" % str(self.targetFeature))
        except AttributeError:
            pass

        return result

    def set_derivedFeature(self, val):
        # string ID to resulting Feature
        self.derivedFeature = val

    def set_targetSeq(self, val):
        # URI to altered sequence record
        self.targetSeq = val

    def set_targetAssembly(self, val):
        # URI to assembly containing altered sequence record
        self.targetAssembly = val

    def set_targetFeature(self, val):
        # string ID to initial Feature
        self.targetFeature = val

    def set_targetRegion(self, val):
        # Region object subject to alteration
        self.targetRegion = val

    def set_agent(self, val):
        # URI for agent that detected this
        self.agent = val

    def set_readCount(self, val):
        self.readCount = val

    def get_type(self):
        return self.type

    def get_derivedFeature(self):
        return self.derivedFeature

    def get_targetSeq(self):
        return self.targetSeq

    def get_targetAssembly(self):
        return self.targetAssembly

    def get_targetFeature(self):
        return self.targetFeature

    def get_targetRegion(self):
        return self.targetRegion

    def get_agent(self):
        return self.agent

    def get_readCount(self):
        return self.readCount



def dumpJSON(obj, fh):
    json.dump(obj, fh, default=jsonobjproc, sort_keys=True, indent=4)


def jsonobjproc(x):
    try:
        return x.__dict__
    except AttributeError:
        try:
            return list(x)
        except TypeError:
            print(f'Unserializable object {x} of type {type(x)}')
            return x


def eng01toenum(eng):
    if eng == 0:
        return ENGINEERING_VALS[0]
    elif eng == 1:
        return ENGINEERING_VALS[2]
    return ENGINEERING_VALS[1]
