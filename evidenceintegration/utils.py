#!/usr/bin/env python3
import logging
import os
import json

class EvidenceIntegrationUtils:

    @staticmethod
    def get_json_dict_from_resources(filename):
        try:
            return json.load((os.path.dirname(os.path.abspath(__file__))
                         ) + "/../examples/" + filename)
        except FileNotFoundError as fnfe:
            logging.error("%s", fnfe)
            pass

    @staticmethod
    def get_json_files_from_examples(sample_id = None):
        jsons = []
        try:
            directory = os.path.dirname(os.path.abspath(__file__)) + "/../resources/"
            for dirpath,_,filenames in os.walk(directory):
                for f in filenames:
                    if f.endswith(".json"):
                        if sample_id is None or sample_id in f:
                            jsons.append(os.path.abspath(os.path.join(dirpath, f)))
        except FileNotFoundError as fnfe:
            logging.debug("%s", fnfe)
        logging.debug("get_json_files_from_resources is returning jsons: %s", str(jsons))
        return jsons


    @staticmethod
    def get_json_files_from_dir(directory, sample_id = None):
        jsons = []
        try:
            for dirpath,_,filenames in os.walk(directory):
                for f in filenames:
                    if f.endswith(".json"):
                        if sample_id is None or sample_id in f:
                            jsons.append(os.path.abspath(os.path.join(dirpath, f)))
        except FileNotFoundError as fnfe:
            logging.debug("%s", fnfe)
        logging.debug("get_json_files_from_dir is returning jsons: %s", str(jsons))
        return jsons