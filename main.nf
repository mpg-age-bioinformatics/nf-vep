#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process get_images {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
    """

    if [[ "${params.containers}" == "singularity" ]] ;

      then

        cd ${params.image_folder}

        if [[ ! -f vep-110.1.sif ]] ;
          then
            singularity pull vep-110.1.sif docker://ensemblorg/ensembl-vep:release_110.1
        fi

        if [[ ! -f rnaseq.python-3.8-2.sif ]] ;
          then
            singularity pull rnaseq.python-3.8-2.sif docker://index.docker.io/mpgagebioinformatics/rnaseq.python:3.8-2
        fi

    fi

    if [[ "${params.containers}" == "docker" ]] ;

      then

        docker pull ensemblorg/ensembl-vep:release_110.1
        docker pull mpgagebioinformatics/rnaseq.python:3.8-2

    fi
    """

}


process get_cache {
  stageInMode 'symlink'
  stageOutMode 'move'

  when:
    ( ! file("${params.genomes}/${params.organism}/${params.release}/vep_cache/${params.organism}/${params.release}_${params.genome_assembly_vep}/info.txt").exists() )

  script:
    """
    INSTALL.pl -c ${params.genomes}/${params.organism}/${params.release}/vep_cache -a cf -s ${params.organism} -y ${params.genome_assembly_vep}
    """
}


process vep {
  stageInMode 'symlink'
  stageOutMode 'move'
  
  input:
    tuple val(pair_id), path(vcf)

  output:
    val pair_id

  when:
    ( ! file("${params.vep_raw_data}/${pair_id}.protein.vcf").exists() ) 

  
  script:
    """
    echo \$(pwd)
    vep --cache -i /workdir/filter/${pair_id}.SNPs.nowt.vcf --species ${params.organism} --biotype --symbol --nearest transcript \
    --dir_cache ${params.genomes}/${params.organism}/${params.release}/vep_cache/ --stats_file /workdir/filter/${pair_id}.protein -o /workdir/filter/${pair_id}.all_variants.vcf

    filter_vep -i /workdir/filter/${pair_id}.all_variants.vcf --filter "BIOTYPE is protein_coding" | filter_vep --filter "IMPACT is MODERATE" |  filter_vep --filter "IMPACT is HIGH" > /workdir/filter/${pair_id}.protein.vcf
    """
}

workflow images {
  main:
    get_images()
}

workflow cache {
  main:
    get_cache()
}

workflow {
  main:
    data = channel.fromFilePairs( "${params.vep_raw_data}/*.SNPs.nowt.vcf", size: -1 )
    vep( data )
}