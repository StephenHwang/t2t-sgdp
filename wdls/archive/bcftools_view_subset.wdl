version 1.0

workflow subset_vcf_sample_list {
    input {
        File inputVCF
        File sample_list
        String population
    }

    call bcftools_view_subset {
        input:
            inputVCF = inputVCF,
            sample_list = sample_list,
            population = population
    }

    call bgzip_bcftools_index {
        input:
            inputVCF = bcftools_view_subset.subsetVCF
    }

    call bcftools_stats {
        input:
            inputVCFgz = bgzip_bcftools_index.bgzipVCF
    }

    output {
        File subsetVCF = bcftools_view_subset.subsetVCF
        File subset_bgzip = bgzip_bcftools_index.bgzipVCF
        File subset_index = bgzip_bcftools_index.index
        File subset_stats = bcftools_stats.stats
    }
}

task bcftools_view_subset {
    input {
        File inputVCF
        File sample_list
        String population
    }

    String vcfPrefix = '~{basename(inputVCF,".vcf")}'

    command <<<
        bcftools view \
            -S "~{sample_list}" \
            --force-samples \
            "~{inputVCF}" > "~{vcfPrefix}.~{population}.vcf"
    >>>

    Int diskGb = ceil(2.0 * size(inputVCF, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ${diskGb} SSD"
        memory: "4G"
        cpu : 2
        preemptible: 3
        maxRetries: 3
    }

    output {
        File subsetVCF = "~{vcfPrefix}.~{population}.vcf"
    }
}

task bgzip_bcftools_index {
    input {
        File inputVCF
    }

    String vcfName = '~{basename(inputVCF)}'

    command <<<
        bgzip -@ "$(nproc)" -c "~{inputVCF}" > "~{vcfName}.gz"
        bcftools index --threads "$(nproc)" "~{vcfName}.gz" -o "~{vcfName}.gz.csi"
    >>>

    Int diskGb = ceil(2.0 * size(inputVCF, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ${diskGb} SSD"
        memory: "4G"
        cpu : 2
        preemptible: 3
        maxRetries: 3
    }

    output {
        File bgzipVCF = "~{vcfName}.gz"
        File index = "~{vcfName}.gz.csi"
    }
}

task bcftools_stats {
    input {
        File inputVCFgz
    }

    String vcfName = '~{basename(inputVCFgz,".vcf.gz")}'

    command <<<
        bcftools stats "~{inputVCFgz}" > "~{vcfName}.bcftools.stats.txt"
    >>>

    Int diskGb = ceil(2.0 * size(inputVCFgz, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ${diskGb} SSD"
        memory: "4G"
        cpu : 2
        preemptible: 3
        maxRetries: 3
    }

    output {
        File stats = "~{vcfName}.bcftools.stats.txt"
    }
}