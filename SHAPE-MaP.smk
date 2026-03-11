configfile : "config.yaml"

rule all:
    input:
        expand("{output_dir}/{sample}/{sample}_{target_name}.shape", sample=config["samples"], target_name=config["target_name"], output_dir=config["output_dir"]),
        expand("{output_dir}/{sample}/{sample}_{target_name}_profiles.pdf", sample=config["samples"], target_name=config["target_name"], output_dir=config["output_dir"]),
        expand("{output_dir}/{sample}/structure/{sample}_{target_name}.dbn", sample=config["samples"], target_name=config["target_name"], output_dir=config["output_dir"]),
        expand("{output_dir}/{sample}/structure/{sample}_{target_name}.ct", sample=config["samples"], target_name=config["target_name"], output_dir=config["output_dir"]),
        expand("{output_dir}/{sample}/structure/{sample}_{target_name}_folding.svg", sample=config["samples"], target_name=config["target_name"], output_dir=config["output_dir"]),


def get_shapemapper_params(wildcards):
    params = [
        "--target", config["target"],
        "--name", wildcards.sample,
        "--out", "{output_dir}/{sample}".format(sample=wildcards.sample, output_dir=config["output_dir"]),
        "--log", "{output_dir}/{sample}/{sample}_shapemapper_log.txt".format(sample=wildcards.sample, output_dir=config["output_dir"]),
        "--verbose",
        "--modified", "--R1", config["samples"][wildcards.sample]["modified"]["R1"], "--R2", config["samples"][wildcards.sample]["modified"]["R2"],
        "--untreated", "--R1", config["samples"][wildcards.sample]["untreated"]["R1"], "--R2", config["samples"][wildcards.sample]["untreated"]["R2"],
        "--output-parsed-mutations",
        "--output-counted-mutations",
        "--per-read-histograms",
    ]
    
    if config["denatured"]:
        params.extend(["--denatured", "--R1", config["samples"][wildcards.sample]["denatured"]["R1"], "--R2", config["samples"][wildcards.sample]["denatured"]["R2"]])
    if config["overwrite"]:
        params.extend(["--overwrite"])
    if config.get("min_depth"):
        params.extend(["--min-depth", str(config["min_depth"])])
    if config.get("max_bg"):
        params.extend(["--max-bg", str(config["max_bg"])])
    if config.get("min_qual_to_trim"):
        params.extend(["--min-qual-to-trim", str(config["min_qual_to_trim"])])
    if config.get("window_to_trim"):
        params.extend(["--window-to-trim", str(config["window_to_trim"])])
    if config.get("min_length_to_trim"):
        params.extend(["--min-length-to-trim", str(config["min_length_to_trim"])])
    if config.get("min_qual_to_count"):
        params.extend(["--min-qual-to-count", str(config["min_qual_to_count"])])
    if config["indiv_norm"]:
        params.extend(["--indiv-norm"])
    
    if config["amplicon"]:
        params.extend(["--amplicon"])
    if config["certain_primer"]:
        params.extend(["--primers", str(config["primers_file"])])
    if config["random_primer"] and config["amplicon"] == False and config["certain_primer"] == False:
        params.extend(["--random-primer-len", str(config["random_primer_length"])])
    if config["amplicon"] or config["certain_primer"]:
        if config.get("max_primer_offset"):
            params.extend(["--max-primer-offset", str(config["max_primer_offset"])])
    return " ".join(params)

rule run_shapemapper:
    input:
        target = config["target"],
        modified_R1 = lambda wildcards: config["samples"][wildcards.sample]["modified"]["R1"],
        modified_R2 = lambda wildcards: config["samples"][wildcards.sample]["modified"]["R2"],
        untreated_R1 = lambda wildcards: config["samples"][wildcards.sample]["untreated"]["R1"],
        untreated_R2 = lambda wildcards: config["samples"][wildcards.sample]["untreated"]["R2"],
    output:
        shape = "{output_dir}/{sample}/{sample}_{target_name}.shape",
        profile = "{output_dir}/{sample}/{sample}_{target_name}_profiles.pdf",
    params:
        cmd_params = get_shapemapper_params,
    threads: 8
    container:
        config["container"]
    shell:
        "shapemapper {params.cmd_params}"

def get_rnastructure_params(wildcards): 
    params = []
    if config.get("temperature"):
        params.extend(["--temperature", str(config["temperature"])])
    if config.get("maximum"):
        params.extend(["--maximum", str(config["maximum"])])
    if config.get("loop"):
        params.extend(["--loop", str(config["loop"])])
    return " ".join(params)


rule predict_structure:
    input:
        shape = "{output_dir}/{sample}/{sample}_{target_name}.shape",
        sequence = config["target"],
    output:
        ct = "{output_dir}/{sample}/structure/{sample}_{target_name}.ct",
    params:
        program_params = get_rnastructure_params
    container:
        config["container"]
    shell:
        "Fold --SHAPE {input.shape} {input.sequence} {output.ct} {params.program_params}"

rule predict_structure_dbn:
    input:
        shape = "{output_dir}/{sample}/{sample}_{target_name}.shape",
        sequence = config["target"],
    output:
        dbn = "{output_dir}/{sample}/structure/{sample}_{target_name}.dbn",
    params:
        program_params = get_rnastructure_params
    container:
        config["container"]
    shell:
        "Fold --SHAPE {input.shape} -k {input.sequence} {output.dbn} {params.program_params}"

rule visualize_structure:
    input:
        shape = "{output_dir}/{sample}/{sample}_{target_name}.shape",
        ct = "{output_dir}/{sample}/structure/{sample}_{target_name}.ct",
    output:
        svg = "{output_dir}/{sample}/structure/{sample}_{target_name}_folding.svg",
    container:
        config["container"]
    shell:
        "draw --SHAPE {input.shape} --svg -n 1 {input.ct} {output.svg}"