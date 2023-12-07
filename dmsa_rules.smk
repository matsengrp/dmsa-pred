build_dir = config.get("build_dir", "builds")

rule prepare_sequences:
    output:
        sequences="data/{lineage}/raw_{segment}.fasta"
    params:
        local_path=config['sequences']
    conda: "../../../workflow/envs/nextstrain.yaml"
    shell:
        """
        cat {params.local_path} | xz -c -d > {output.sequences}
        """

# TODO how does metadata work?
# looks like select strains is the smk for parsing fasta fields.
rule prepare_metadata:
    output:
        # titers="data/{lineage}/raw_{segment}.tsv"
        titers="data/{lineage}/metadata_{segment}.tsv"
    params:
        local_path=config['metadata']
    conda: "../../../workflow/envs/nextstrain.yaml"
    shell:
        """
        cat {params.local_path} | xz -c -d > {output.titers}
        """

import os
# predict escape on the HA segment translations.
rule variant_escape_prediction:
    input:
        translations=build_dir + "/{build_name}/ha/translations.done",
    output:
        node_data = build_dir + "/{build_name}/{segment}/dmsa-phenotype/{collection}/{experiment}_variant_escape_prediction.json",
        pred_data = build_dir + "/{build_name}/{segment}/dmsa-phenotype/{collection}/{experiment}_variant_escape_prediction.csv",
    log:
        "logs/{build_name}/{segment}/{collection}/{experiment}_variant_escape_prediction.txt"
    params:
        basedir = lambda w: config["dmsa_phenotype_collections"].get(w.collection)['mut_effects_dir'],
        dms_wt_seq_id = lambda w: config["dmsa_phenotype_collections"].get(w.collection)['dms_wt_seq_id'],
        mut_effect_col = lambda w: config["dmsa_phenotype_collections"].get(w.collection)['mut_effect_col'],
        mutation_col = lambda w: config["dmsa_phenotype_collections"].get(w.collection)['mutation_col'],
        allow_unmeasured_aa_subs_at_these_sites = lambda w: config["dmsa_phenotype_collections"].get(w.collection)['allow_unmeasured_aa_subs_at_these_sites'],
        mut_effects_df = lambda w: os.path.join(
            config["dmsa_phenotype_collections"].get(w.collection)['mut_effects_dir'], 
            w.experiment
        ),
        alignment = lambda w: build_dir + f"/{w.build_name}/{w.segment}/nextalign/masked.gene.HA_withInternalNodes.fasta"
    conda:
        "dmsa_env.yaml"
    resources:
        mem_mb=2000
    shell:
        """
        python profiles/dmsa-phenotype/dmsa-pred/dmsa_pred.py phenotype-prediction \
            --model-type additive \
            --alignment {params.alignment} \
            --dms-wt-seq-id {params.dms_wt_seq_id} \
            --mut-effects-df {params.mut_effects_df} \
            --mut-effect-col {params.mut_effect_col} \
            --mutation-col {params.mutation_col} \
            --allow-unmeasured-aa-subs-at-these-sites {params.allow_unmeasured_aa_subs_at_these_sites} \
            --experiment-label {wildcards.experiment} \
            --output-json {output.node_data} \
            --output-df {output.pred_data} 2>&1 | tee {log}
        """

# import glob
# # define the required inputs for the export rule
# def _get_escape_node_data_by_wildcards(wildcards):
#     """Return a list of node data files to include for a given build's wildcards.
#     """
#     # Define inputs shared by all builds.
#     wildcards_dict = dict(wildcards)

#     # TODO you could check that the 'HA' segment is in SEGMENTS['ha']?
#     if not config.get("dmsa_phenotype_collections", False):
#         return []

#     build = config["builds"][wildcards.build_name]
#     if build.get('dmsa_phenotype', False) and wildcards.segment == 'ha':
#         for collection in build.get('dmsa_phenotype'):

#             # get the kwargs for the collection of dms escape models
#             kwargs = config["dmsa_phenotype_collections"][collection]

#             # run the predictions using every csv in the glob
#             requested_files = expand(
#                 rules.variant_escape_prediction.output.node_data,
#                 collection=collection,
#                 experiment=[
#                     os.path.basename(fp) 
#                     for fp in glob.glob(kwargs['mut_effects_dir']+"/*.csv")
#                 ],
#                 **wildcards_dict
#             )
#             print("\n".join(requested_files))
#             inputs.extend(requested_files)

#     return inputs

# # # overwrite the export rule to include escape node data
# rule export:
#     message: "Exporting data files for auspice"
#     input:
#         tree = rules.refine.output.tree,
#         metadata = build_dir + "/{build_name}/metadata.tsv",
#         node_data = _get_node_data_by_wildcards + _get_escape_node_data_by_wildcards,
#         auspice_config = lambda w: config['builds'][w.build_name]['auspice_config'],
#         lat_longs = config['lat-longs']
#     output:
#         auspice_json = "auspice/{build_name}_{segment}.json",
#         root_sequence_json = "auspice/{build_name}_{segment}_root-sequence.json",
#     conda: "../../envs/nextstrain.yaml"
#     benchmark:
#         "benchmarks/export_{build_name}_{segment}.txt"
#     log:
#         "logs/export_{build_name}_{segment}.txt"
#     shell:
#         """
#         augur export v2 \
#             --tree {input.tree} \
#             --metadata {input.metadata} \
#             --node-data {input.node_data} \
#             --include-root-sequence \
#             --lat-longs {input.lat_longs} \
#             --auspice-config {input.auspice_config} \
#             --output {output.auspice_json} 2>&1 | tee {log}
#         """
