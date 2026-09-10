#!/home/liuzm/miniconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 19:01:05 2022

Author: Zhang Zhichao
Email: zhichao.zhang@glbizzia.com

"""

import os
import argparse
from argparse import RawTextHelpFormatter


bin_path = "/gpfs/users/liuzm/software/SingleCell/data"
go_base = "%s/go_base" % bin_path
ko2map_file = "%s/k_id2map" % bin_path
kegg_pathway_file = "%s/kegg_pathway" % bin_path
kegg_all_organs_file = "%s/all_kegg_organs" % bin_path
species_pathway_path = "%s/kegg_data" % bin_path


parser = argparse.ArgumentParser(prog="eggnog_anno.py", description="""eggnogmapper annotation pipline ver 1.0

""", formatter_class=RawTextHelpFormatter)
parser.add_argument("-t", "--type", dest="type", required=True, choices=["inhouse", "eggnog5"], help="choice input file type.")
parser.add_argument("-g", "--go", dest="go", required=False, help="go annotation file, require if type was inhouse.")
parser.add_argument("-k", "--kegg", dest="kegg", required=False, help="KEGG annotation file, require if type was inhouse.")
parser.add_argument("-s", "--species", dest="species", default=None, required=False, help="KEGG species abbr, %s " % kegg_all_organs_file)
parser.add_argument("-e", "--eggnog5", dest="eggnog5", required=False, help="eggnog5 annotation file, require if eggnog5 type was chosen.")
# parser.add_argument("-t", "--type", dest="type", required=True, choices=["inhouse", "eggnog", "custom"], help="choice input file type.")
parser.add_argument("-o", "--output", dest="output", required=False, help="shell out put path, default ./ ")


args = parser.parse_args()
absp = os.path.abspath
opj = os.path.join

# in_pepfile = absp(args.pep)
in_go_file = args.go
in_kegg_file = args.kegg
species = args.species
in_eggnog_file = args.eggnog5
pipline_type = args.type
outdir = args.output


if outdir == None:
    outdir = os.getcwd()

go_output_file_path = opj(outdir, "GO.universal.tsv")
kegg_output_file_path = opj(outdir, "KEGG.universal.tsv")


go_output_file_title = "Query\tGO\tDescription\tLevel\n"
kegg_output_file_title = "Query\tKO\tPathway\tDescription\n"


if species == None:
    species_prefix = "map"
else:
    species_prefix = species


def parse_go_database():
    go_base_dict = {}
    with open(go_base) as go_database:
        for line in go_database:
            line = line.strip().split("\t")
            go_term_id = line[0]
            go_discription = line[1]
            go_level = line[2]
            go_base_dict[go_term_id] = f"{go_term_id}\t{go_discription}\t{go_level}"
        return go_base_dict


def spec_info():
    with open(kegg_all_organs_file) as kegg_all_organs:
        for line in kegg_all_organs:
            line = line.strip().split("\t")
            kegg_abbr = line[1]
            if species == kegg_abbr:
                return line
    return None


def species_pathway():
    species_pathway_list = []
    species_pathway_file_path = opj(species_pathway_path, species)

    if os.path.isfile(species_pathway_file_path):
        species_info = spec_info()
        if species_info:
            note = "Species info: %s" % "\t".join(species_info)
            print(note)
        else:
            note = "Not found %s in KEGG organs list."
            print(note)

    else:
        os.system('echo -e "\033[31m Wrong KEGG species abbr. Please check. \033"')
        raise ValueError("Wrong KEGG species abbr. Please check.")

    with open(species_pathway_file_path) as species_pathway_file:
        for line in species_pathway_file:
            species_pathway_list.append(line.strip())
    return species_pathway_list


def parse_kegg_database():
    ko2map_dict = {}
    kegg_pathway_dict = {}

    if species != None:
        species_pathway_list = species_pathway()

    with open(ko2map_file) as in_ko2map, open(kegg_pathway_file) as kegg_pathway:
        for line in in_ko2map:
            line = line.strip().split("\t")
            kegg_k_id = line[0]
            pathway_id = line[1]

            if species != None:
                try:
                    species_pathway_list.index(pathway_id)

                except ValueError:
                    continue

            try:
                ko2map_dict[kegg_k_id].append(pathway_id)

            except KeyError:
                ko2map_dict[kegg_k_id] = [pathway_id]

        for line in kegg_pathway:
            line = line.strip().split("\t")
            pathway_id = line[0]
            pathway_desc = line[1]
            kegg_pathway_dict[pathway_id] = pathway_desc

        return ko2map_dict, kegg_pathway_dict


def parse_inhouse():
    go_base_dict = parse_go_database()
    ko2map_dict, kegg_pathway_dict = parse_kegg_database()
    with open(in_go_file) as go_file, open(in_kegg_file) as kegg_file, open(go_output_file_path, "w") as go_output_file, open(kegg_output_file_path, "w") as kegg_output_file:
        go_output_file.write(go_output_file_title)
        kegg_output_file.write(kegg_output_file_title)

        for line in go_file:
            line = line.strip().split("\t")
            gene_id = line[0]
            go_term_list = line[1:]
            for go_term in go_term_list:
                try:
                    new_line = "%s\t%s\n" % (gene_id, go_base_dict[go_term])
                    go_output_file.write(new_line)

                except:
                    message = "%s was not found in the GO.base" % go_term
                    print(message)

        for line in kegg_file:
            line = line.strip().split("\t")
            gene_id = line[0]
            kegg_k_id_list = line[1].split(",")

            for k_id in kegg_k_id_list:
                try:
                    pathway_list = [[x, kegg_pathway_dict[x]] for x in ko2map_dict[k_id]]
                    for pathway_list_tmp in pathway_list:
                        pathway_id = pathway_list_tmp[0].replace("map", species_prefix)
                        description = pathway_list_tmp[1]
                        new_line = "%s\t%s\t%s\t%s\n" % (gene_id, k_id, pathway_id, description)
                        kegg_output_file.write(new_line)

                except:
                    pass


def parse_eggnog():
    go_base_dict = parse_go_database()
    ko2map_dict, kegg_pathway_dict = parse_kegg_database()

    with open(in_eggnog_file) as in_eggnog, open(go_output_file_path, "w") as go_output_file, open(kegg_output_file_path, "w") as kegg_output_file:
        go_output_file.write(go_output_file_title)
        kegg_output_file.write(kegg_output_file_title)

        for line in in_eggnog:
            if not line.startswith("#"):
                line = line.strip().split("\t")
                gene_id = line[0]
                go_term_list = line[9].split(",")
                if go_term_list == ["-"]:
                    pass
                else:
                    for go_term in go_term_list:
                        try:
                            new_line = "%s\t%s\n" % (gene_id, go_base_dict[go_term])
                            go_output_file.write(new_line)

                        except:
                            message = "%s was not found in the GO.base" % go_term
                            print(message)

                kegg_k_id_list = line[11].replace("ko:", "").split(",")
                if kegg_k_id_list == ["-"]:
                    pass

                else:
                    for k_id in kegg_k_id_list:
                        try:
                            pathway_list = [[x, kegg_pathway_dict[x]] for x in ko2map_dict[k_id]]
                            for pathway_list_tmp in pathway_list:
                                pathway_id = pathway_list_tmp[0].replace("map", species_prefix)
                                description = pathway_list_tmp[1]
                                new_line = "%s\t%s\t%s\t%s\n" %(gene_id, k_id, pathway_id, description)
                                kegg_output_file.write(new_line)

                        except:
                            pass


def main():
    if pipline_type == "inhouse":
        if in_go_file == None or in_kegg_file == None:
            raise ValueError("""Please input both go and kegg file with parameters "--go and --kegg" """)
        parse_inhouse()

    elif pipline_type == "eggnog5":
        if in_eggnog_file == None:
            raise ValueError("""Please input eggnog5 annotation file with parameter "--eggnog5" """)

        parse_eggnog()


    # elif pipline_type == "custom":
    #     parser_type = "custom"
    #     if in_go_file == None or in_kegg_file == None:
    #         raise ValueError("Please input both go and kegg file")


    else:
        raise ValueError("Wrong pipline type, only support inhouse or eggnog pipline function file")



if __name__ == '__main__':
    main()




