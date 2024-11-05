from pathlib import Path
from datetime import datetime
from argparse import ArgumentParser
import json, toml, pandas as pd
from csv import QUOTE_ALL
from rich.console import Console
from rich_argparse import RichHelpFormatter

from typing import List, Dict


class Converter:

    transnetyx_config_filename = "transnetyx_config.json"
    pyrat_config_filename = "pyrat_config.toml"

    file_export_prefix = "PYRAT_EXPORT"

    def __init__(self, file_name: str, *, root_folder=Path.home() / "Downloads"):

        now = datetime.now().strftime(r"%Y_%m_%dT%H-%M-%S")

        self.input_fullpath = Path(root_folder) / file_name
        self.output_fullpath = Path(root_folder) / f"{self.file_export_prefix}_{file_name.split(".")[0]}_{now}.csv"

    def load_config(self):
        self.transnetyx_config = self.get_transnetyx_config()
        self.pyrat_config = self.get_pyrat_config()

    def get_transnetyx_config(self):
        file_path = self.resolve_path(self.transnetyx_config_filename)
        with open(file_path, "r") as file:
            return json.load(file)

    def get_pyrat_config(self):
        file_path = self.resolve_path(self.pyrat_config_filename)
        with open(file_path, "r") as file:
            return toml.load(file)

    def resolve_path(self, path):
        base_folder = Path(__file__).parent
        path = Path(path)
        if path.is_absolute():
            return path
        else:
            return base_folder.joinpath(path).resolve()

    def get_transnetyx_data(self):
        return pd.read_csv(self.input_fullpath)

    def get_pyrat_data(self, transnetyx_data):
        return pd.DataFrame(transnetyx_data.apply(self.get_mutations, axis="columns").tolist())  # .replace(pd.NA, "")

    def write_pyrat_data(self, pyrat_data: pd.DataFrame):
        pyrat_data.to_csv(self.output_fullpath, sep="\t", index_label=False, index=False, quoting=QUOTE_ALL, na_rep="")

    def convert(self):
        self.load_config()
        self.write_pyrat_data(self.get_pyrat_data(self.get_transnetyx_data()))

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
            mutations_columns[f"genotype {index + 1}"] = self.get_genotype(mutation_map, animal_data)

        return dict(id=id, line=pyrat_strain, **mutations_columns)

    def get_genotype_string(self, allele_list: List["AlleleValue"], reverse=True):
        return "/".join([g.to_marker() for g in sorted(allele_list, key=lambda v: bool(v), reverse=reverse)])

    def get_genotype(self, mutation_map: Dict, animal_data: pd.Series):

        allele_list = []
        for allele in mutation_map["alleles"]:
            if isinstance(allele, str):
                allele_list.append(AlleleValue(allele, **mutation_map))
            else:
                presence = animal_data[allele["name"]]
                allele_blueprint = AlleleValue(allele["if_positive"], **mutation_map)
                real_allele = allele_blueprint.check_presence(presence)
                allele_list.append(real_allele)

        return self.get_genotype_string(allele_list, reverse=AlleleValue.is_mutant(mutation_map["leading_marker"]))


class AlleleValue:

    mutant_tag = "mut"
    wild_type_tag = "wt"

    @staticmethod
    def is_mutant(value):
        return True if value == AlleleValue.mutant_tag else False

    def __init__(self, value, *, mut_marker, wt_marker, **_):
        self.value = value
        self.mut_marker = mut_marker
        self.wt_marker = wt_marker

    @property
    def kwargify(self):
        return {"mut_marker": self.mut_marker, "wt_marker": self.wt_marker}

    def __neg__(self):
        return AlleleValue(self.wild_type_tag if self.is_mutant(self.value) else self.mutant_tag, **self.kwargify)

    def __str__(self):
        return self.value

    def __repr__(self):
        return self.__str__()

    def check_presence(self, presence_marker):
        if "+" in presence_marker:
            return self
        return -self

    def to_marker(self):
        return self.mut_marker if self.is_mutant(self.value) else self.wt_marker

    def __bool__(self):
        return self.is_mutant(self.value)


def run():

    console = Console()
    parser = ArgumentParser(formatter_class=RichHelpFormatter)
    parser.add_argument(
        "file",
        nargs="?",
        help="The file name or path leading to the file to be converted. By default, the file name is searched for in "
        "the Downloads folder of the current user's home.",
    )
    parser.add_argument(
        "-f",
        "--file",
        dest="file",
        help="The file name or path leading to the file to be converted. By default, the file name is searched for in "
        "the Downloads folder of the current user's home.",
    )
    args = parser.parse_args()

    if not args.file:
        parser.error("You must supply a filename either as first argument or with the -f or --file flag.")

    Converter(args.file).convert()

    console.print("🐭 Convertion done ! ✅", style="bright_green bold")
