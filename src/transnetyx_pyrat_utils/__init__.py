from pathlib import Path
from datetime import datetime
import json, toml, pandas as pd
from csv import QUOTE_ALL

from typing import List, Dict


class Converter:

    transnetyx_config_filename = "transnetyx_config.toml"
    pyrat_config_filename = "pyrat_config.toml"

    def __init__(self, file_name: str, *, root_folder=Path.home()):

        now = datetime.now().strftime(r"%Y_%m_%dT%H-%M-%S")

        self.input_fullpath = Path(root_folder) / file_name
        self.output_fullpath = f"PYRAT_EXPORT_{file_name.split(".")[0]}_{now}.csv"

        self.load_config()

    def load_config(self):
        self.transnetyx_config = self.get_transnetyx_config()
        self.pyrat_config = self.get_pyrat_config()

    def get_transnetyx_config(self):
        file_path = Path(self.transnetyx_config_filename).resolve()
        with open(file_path, "r") as file:
            return json.load(file)

    def get_pyrat_config(self):
        file_path = Path(self.pyrat_config_filename).resolve()
        with open(file_path, "r") as file:
            return toml.load(file)

    def get_transnetyx_data(self):
        return pd.read_csv(self.input_fullpath)

    def write_pyrat_data(self, pyrat_data: pd.DataFrame):
        pyrat_data.to_csv(self.output_fullpath, sep="\t", index_label=False, index=False, quoting=QUOTE_ALL, na_rep="")

    def convert(self, transnetyx_data):
        return pd.DataFrame(transnetyx_data.apply(self.get_mutations, axis="columns").tolist())  # .replace(pd.NA, "")

    def get_mutations(self, animal_data):

        id_number_zfill = self.pyrat_config["animals"]["id_zfill"]
        id_prefix = self.pyrat_config["animals"]["id_prefix"]
        id_number = animal_data["Sample"].split(" ")[0].zfill(id_number_zfill)
        id = f"{id_prefix}{id_number}"

        transnetyx_strain = animal_data["Strain"]
        pyrat_strain = self.transnetyx_config[transnetyx_strain]["pyrat_name"]
        mutations = self.transnetyx_config[transnetyx_strain]["mutations"]

        mutations_columns = {}
        for index, mutation_map in enumerate(mutations):
            mutations_columns[f"mutation {index + 1}"] = mutation_map["name"]
            mutations_columns[f"genotype {index + 1}"] = get_genotype(mutation_map, animal_data)

        return dict(id=id, line=pyrat_strain, **mutations_columns)


class AlleleValue:

    mutant_tag = "mut"
    wild_type_tag = "wt"

    def __init__(self, value, *, mut_marker, wt_marker, **_):
        self.value = value
        self.mut_marker = mut_marker
        self.wt_marker = wt_marker

    @property
    def kwargify(self):
        return {"mut_marker": self.mut_marker, "wt_marker": self.wt_marker}

    def __neg__(self):
        return AlleleValue(self.mutant_tag if self.value == self.wild_type_tag else self.wild_type_tag, **self.kwargify)

    def __str__(self):
        return self.value

    def __repr__(self):
        return self.__str__()

    def check_presence(self, presence_marker):
        if "+" in presence_marker:
            return self
        return -self

    def to_marker(self):
        return self.mut_marker if self.value == self.mutant_tag else self.wt_marker

    def __bool__(self):
        return True if self.value == self.mutant_tag else False


def get_genotype_string(allele_list: List[AlleleValue]):
    return "/".join([g.to_marker() for g in sorted(allele_list, key=lambda v: bool(v), reverse=True)])


def get_genotype(mutation_map: Dict, animal_data: pd.Series):

    allele_list = []
    for allele in mutation_map["alleles"]:
        if isinstance(allele, str):
            allele_list.append(AlleleValue(allele, **mutation_map))
        else:
            presence = animal_data[allele["name"]]
            allele_blueprint = AlleleValue(allele["if_positive"], **mutation_map)
            real_allele = allele_blueprint.check_presence(presence)
            allele_list.append(real_allele)

    return get_genotype_string(allele_list)
