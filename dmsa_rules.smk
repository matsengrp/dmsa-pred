rule variant_escape_fraction_prediction:
    input:
        alignments = "results/{build_name}/translations/aligned.gene.S_withInternalNodes.fasta",
    output:
        node_data = "results/{build_name}/escape_fraction_pred.json"
    log:
        "logs/variant_escape_fraction_prediction_{build_name}.txt"
    params:
        dms_wt_seq = lambda w: config["mutation_escape_fractions"]["wt_seq"],
        mut_effects_df = lambda w: config["mutation_escape_fractions"]["mut_effects_df"],
    conda:
        config["conda_environment"],
    resources:
        mem_mb=2000
    shell:
        """
        python3 scripts/escape_frac_predict.py \
            --alignment {input.alignments} \
            --mut-effects-df {params.mut_effects_df} \
            --dms-wt-seq {params.dms_wt_seq} \
            --output {output.node_data} 2>&1 | tee {log}
        """

rule polyclonal_escape_prediction:
    input:
        alignments = "results/{build_name}/translations/aligned.gene.S_withInternalNodes.fasta",
    output:
        node_data = "results/{build_name}/{antibody}/escape_pred.json"
    log:
        "logs/polyclonal_escape_prediction_{build_name}_{antibody}.txt"
    params:
        dms_wt_seq = lambda w: config["polyclonal_antibody_models"][f"{w.antibody}"]["wt_seq"],
        activity_wt_df = lambda w: config["polyclonal_antibody_models"][f"{w.antibody}"]["activity_wt_df"],
        mut_escape_df = lambda w: config["polyclonal_antibody_models"][f"{w.antibody}"]["mut_escape_df"],
        concentrations = lambda w: config["polyclonal_antibody_models"][f"{w.antibody}"]["concentrations"]
    conda:
        "dmsa_env.yaml"
    resources:
        mem_mb=2000
    shell:
        """
        dmsa_pred/dmsa_pred.py polyclonal-escape \
            --alignment {input.alignments} \
            --mut-escape-df {params.mut_escape_df} \
            --activity-wt-df {params.activity_wt_df} \
            --dms-wt-seq-id {params.dms_wt_seq} \
            --serum-label {wildcards.antibody} \
            --concentrations-list {params.concentrations} \
            --output {output.node_data} 2>&1 | tee {log}
        """
