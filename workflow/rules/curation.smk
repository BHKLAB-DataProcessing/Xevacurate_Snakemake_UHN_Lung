rule scan:
    output:
        json="results/scan/all_file_scan.json"
    params:
        raw_dir="data/input/drug_screen",
        out_dir="results/scan",
        remove_list="data/input/additional/models_to_remove.xlsx"  # change/keep; script can ignore if missing
    log:
        "logs/scan/scan.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.out_dir}" "logs/scan"
        PYTHONPATH=src pixi run python workflow/scripts/run_scan.py \
          --drug-screen-dir "{params.raw_dir}" \
          --output-dir "{params.out_dir}" \
          --remove-list "{params.remove_list}" \
          > "{log}" 2>&1
        """

rule process:
    input:
        manifest="results/scan/all_file_scan.json"
    output:
        tsv="results/process/processed_data.tsv",
        qc="results/process/qc_report.json"
    params:
        base_dir="data/input/drug_screen",
        out_dir="results/process"
    log:
        "logs/process/process.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.out_dir}" "logs/process"
        PYTHONPATH=src pixi run python workflow/scripts/run_process.py \
          --manifest "{input.manifest}" \
          --base-dir "{params.base_dir}" \
          --output-dir "{params.out_dir}" \
          > "{log}" 2>&1
        """

rule update:
    input:
        tsv="results/process/processed_data.tsv"
    output:
        batch="results/update/update_model_drug_batch_data.tsv",
        model="results/update/update_model_data.tsv",
        drug="results/update/update_model_drug_data.tsv",
        filtered="results/update/drug_filtered_out.tsv"
    params:
        out_dir="results/update",
        model_map="data/input/additional/all_model_clinical.csv",
        drug_map="data/input/additional/all_drug_mapping.xlsx"
    log:
        "logs/update/update.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.out_dir}" "logs/update"
        PYTHONPATH=src pixi run python workflow/scripts/run_update.py \
          --processed-tsv "{input.tsv}" \
          --model-map-file "{params.model_map}" \
          --drug-map-file "{params.drug_map}" \
          --out-dir "{params.out_dir}" \
          > "{log}" 2>&1
        """

rule build:
    input:
        batch="results/update/update_model_drug_batch_data.tsv"
    output:
        "results/build/model.csv",
        "results/build/drug.csv",
        "results/build/experiment.csv",
        "results/build/expDesign.csv",
        "results/build/modToBiobaseMap.csv",
        "results/build/CNV_pdata.csv",
        "results/build/mutation_pdata.csv",
        "results/build/RNASeq_pdata.csv"
    params:
        out_dir="results/build",
        omics_dir="data/input/omics"
    log:
        "logs/build/build.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.out_dir}" "logs/build"
        PYTHONPATH=src pixi run python workflow/scripts/run_build.py \
          --in "{input.batch}" \
          --out-dir "{params.out_dir}" \
          --omics-dir "{params.omics_dir}" \
          > "{log}" 2>&1
        """

rule create_xevaset:
    input:
        expand(
            "results/build/{f}.csv",
            f=[
                "model", "drug", "experiment", "expDesign",
                "modToBiobaseMap", "CNV_pdata", "mutation_pdata", "RNASeq_pdata"
            ],
        )
    output:
        # derive from config to avoid duplication
        rds=lambda wildcards: "results/xevaset/{}".format(
            config.get("xevaset", {}).get(
                "rds_name", "UHN_Tsao_Lung_DrugResponse_2022_v1.rds"
            )
        )
    params:
        image     = config.get("docker", {}).get("xeva_r_image", "docker.io/gfeng2023/xeva-build-r:bioc-3.20"),
        build_dir = "results/build",
        omics_dir = "data/input/omics",
        out_dir   = "results/xevaset",
        name      = config.get("xevaset", {}).get("name", "Tsao_Lung_2022"),
        rds_name  = config.get("xevaset", {}).get("rds_name", "UHN_Tsao_Lung_DrugResponse_2022_v1.rds"),
        # uncomment if you want to force a specific arch (e.g., on Apple Silicon running amd64 layers):
        # platform = "--platform linux/amd64"
        platform = ""
    log:
        "logs/xevaset/xevaset.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.out_dir}" "logs/xevaset"

        # run inside container, write files as host user, capture ALL output to log
        {{
          docker run --rm \
            -u "$(id -u):$(id -g)" \
            -v "$PWD":/work -w /work \
            {params.platform} \
            {params.image} \
            Rscript workflow/scripts/create_xevaset.R \
              --build-dir {params.build_dir} \
              --omics-dir {params.omics_dir} \
              --out-dir {params.out_dir} \
              --name {params.name} \
              --rds-name {params.rds_name}
        }} > "{log}" 2>&1

        # sanity check
        test -s "{output.rds}"
        """
