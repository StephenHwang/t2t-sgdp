version 1.0

workflow haplotype_calling_chrom {
	input {
		File refFasta
		File fastaIndex
		File fastaDict
		File cram
		File cramIndex
		File? parXbed
		File? nonparXbed
		File? parYbed
		String sampleName
		String chrXname = "chrX"
		String chrYname = "chrY"
		String sex
	}

	if (sex == "female") {
		call diploidHC as hc_XX_X {
			input:
				refFasta = refFasta,
				fastaIndex = fastaIndex,
				fastaDict = fastaDict,
				cram = cram,
				cramIndex = cramIndex,
				chrom = chrXname,
				sampleName = sampleName
		}
	}

	if (sex == "male") {
		call sexchrHC as hc_XY_X_nonPAR {
			input:
				refFasta = refFasta,
				fastaIndex = fastaIndex,
				fastaDict = fastaDict,
				cram = cram,
				cramIndex = cramIndex,
				region = chrXname,
				regionName = "chrX_non_PAR",
				excludeInterval = parXbed,
				ploidy = 1,
				sampleName = sampleName
		}
		call sexchrHC as hc_XY_X_PAR {
			input:
				refFasta = refFasta,
				fastaIndex = fastaIndex,
				fastaDict = fastaDict,
				cram = cram,
				cramIndex = cramIndex,
				region = chrXname,
				regionName = "chrX_PAR",
				excludeInterval = nonparXbed,
				ploidy = 2,
				sampleName = sampleName
		}
		call sexchrHC as hc_XY_Y_nonPAR {
			input:
				refFasta = refFasta,
				fastaIndex = fastaIndex,
				fastaDict = fastaDict,
				cram = cram,
				cramIndex = cramIndex,
				region = chrYname,
				regionName = "chrY",
				excludeInterval = parYbed,
				ploidy = 1,
				sampleName = sampleName
		}
	}

	output {
		File? XX_X_hcVCF = hc_XX_X.hcVCF
		File? XX_X_hcVCF_gz = hc_XX_X.hcVCF_gz
		File? XX_X_hcVCF_gz_tbi = hc_XX_X.hcVCF_gz_tbi
		File? XY_X_non_PAR_hcVCF = hc_XY_X_nonPAR.hcVCF
		File? XY_X_non_PAR_hcVCF_gz = hc_XY_X_nonPAR.hcVCF_gz
		File? XY_X_non_PAR_hcVCF_gz_tbi = hc_XY_X_nonPAR.hcVCF_gz_tbi
		File? XY_X_PAR_hcVCF = hc_XY_X_PAR.hcVCF
		File? XY_X_PAR_hcVCF_gz = hc_XY_X_PAR.hcVCF_gz
		File? XY_X_PAR_hcVCF_gz_tbi = hc_XY_X_PAR.hcVCF_gz_tbi
		File? XY_Y_nonPAR_hcVCF = hc_XY_Y_nonPAR.hcVCF
		File? XY_Y_nonPAR_hcVCF_gz = hc_XY_Y_nonPAR.hcVCF_gz
		File? XY_Y_nonPAR_hcVCF_gz_tbi = hc_XY_Y_nonPAR.hcVCF_gz_tbi
	}
}

task diploidHC {
	input {
		File refFasta
		File fastaIndex
		File fastaDict
		File cram
		File cramIndex
		String chrom
		String sampleName
	}

	String cramName = '~{basename(cram)}'
	String fastaName = '~{basename(refFasta)}'
	String refBaseFasta = '~{basename(refFasta, ".fasta")}'
	String refBase = '~{basename(refBaseFasta, ".fa")}'

	command <<<
		mkdir -p inputs/
		mv "~{cram}" inputs/
		mv "~{cramIndex}" inputs/

		mv "~{refFasta}" inputs/
		mv "~{fastaIndex}" inputs/
		mv "~{fastaDict}" "inputs/~{refBase}.dict"

		gatk HaplotypeCaller \
			--java-options "-Xmx32G -XX:+UseParallelGC -XX:ParallelGCThreads=$(nproc) -Djava.io.tmpdir=/dev/shm" \
			-R inputs/"~{fastaName}" \
			-I inputs/"~{cramName}" \
			-L "~{chrom}" \
			-pairHMM AVX_LOGLESS_CACHING \
			-O "~{sampleName}.~{chrom}.hc.vcf" \
			-ERC GVCF \
			-ploidy 2 \
			-A Coverage \
			-A DepthPerAlleleBySample \
			-A DepthPerSampleHC \
			-A InbreedingCoeff \
			-A MappingQualityRankSumTest \
			-A MappingQualityZero \
			-A QualByDepth \
			-A ReadPosRankSumTest \
			-A RMSMappingQuality \
			-A StrandBiasBySample

		bgzip -c "~{sampleName}.~{chrom}.hc.vcf" > "~{sampleName}.~{chrom}.hc.vcf.gz"
		tabix "~{sampleName}.~{chrom}.hc.vcf.gz"
	>>>

	Int diskGb = ceil(2.0 * size(cram, "G"))

	runtime {
		docker : "szarate/t2t_variants:v0.0.2"
		disks : "local-disk ${diskGb} SSD"
		memory: "24G"
		cpu : 2
		preemptible: 3
		maxRetries: 3
	}

	output {
		File hcVCF = "~{sampleName}.~{chrom}.hc.vcf"
		File hcVCF_gz = "~{sampleName}.~{chrom}.hc.vcf.gz"
		File hcVCF_gz_tbi = "~{sampleName}.~{chrom}.hc.vcf.gz.tbi"
	}
}

task sexchrHC {
	input {
		File refFasta
		File fastaIndex
		File fastaDict
		File cram
		File cramIndex
		String region
		File? regionBed
		String regionName
		File? excludeInterval
		String sampleName
		Int ploidy
	}

	String excludeName = '~{(if defined(excludeInterval) then basename(select_first([excludeInterval])) else "")}'
	String regions = '~{(if defined(regionBed) then basename(select_first([regionBed])) else "")}'

	String cramName = '~{basename(cram)}'
	String fastaName = '~{basename(refFasta)}'
	String refBaseFasta = '~{basename(refFasta, ".fasta")}'
	String refBase = '~{basename(refBaseFasta, ".fa")}'

	command <<<
		mkdir -p inputs/

		mv "~{cram}" inputs/
		mv "~{cramIndex}" inputs/

		mv "~{refFasta}" inputs/
		mv "~{fastaIndex}" inputs/
		mv "~{fastaDict}" "inputs/~{refBase}.dict"

		echo "~{excludeName}"
		echo "~{regions}"

		if [[ -n "~{regions}" ]]; then
			mv "~{regionBed}" inputs/
			specifyRegion="--intervals inputs/~{regions}"
		else
			specifyRegion="--intervals ~{region}"
		fi

		if [[ -n "~{excludeName}" ]]; then
			mv "~{excludeInterval}" inputs/
			specifyExclude="--exclude-intervals inputs/~{excludeName}"
		else
			specifyExclude=""
		fi

		gatk HaplotypeCaller \
			--java-options "-Xmx32G -XX:+UseParallelGC -XX:ParallelGCThreads=$(nproc) -Djava.io.tmpdir=/dev/shm" \
			-R inputs/"~{fastaName}" \
			-I inputs/"~{cramName}" \
			${specifyRegion} ${specifyExclude} \
			-pairHMM AVX_LOGLESS_CACHING \
			-O "~{sampleName}.~{regionName}.hc.vcf" \
			-ERC GVCF \
			-ploidy "~{ploidy}" \
			-A Coverage \
			-A DepthPerAlleleBySample \
			-A DepthPerSampleHC \
			-A InbreedingCoeff \
			-A MappingQualityRankSumTest \
			-A MappingQualityZero \
			-A QualByDepth \
			-A ReadPosRankSumTest \
			-A RMSMappingQuality \
			-A StrandBiasBySample
		
		bgzip -c "~{sampleName}.~{regionName}.hc.vcf" > "~{sampleName}.~{regionName}.hc.vcf.gz"
		tabix "~{sampleName}.~{regionName}.hc.vcf.gz"
	>>>

	Int diskGb = ceil(2.0 * size(cram, "G"))

	runtime {
		docker : "szarate/t2t_variants:v0.0.2"
		disks : "local-disk ${diskGb} SSD"
		memory: "24G"
		cpu : 2
		preemptible: 3
		maxRetries: 3
	}

	output {
		File hcVCF = "~{sampleName}.~{regionName}.hc.vcf"
		File hcVCF_gz = "~{sampleName}.~{regionName}.hc.vcf.gz"
		File hcVCF_gz_tbi = "~{sampleName}.~{regionName}.hc.vcf.gz.tbi"
	}
}
