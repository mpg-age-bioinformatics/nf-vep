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
    vep --cache -i /workdir/filter/${pair_id}.SNPs.nowt.vcf --species ${params.organism} --biotype --symbol --nearest transcript --vcf \
    --dir_cache ${params.genomes}/${params.organism}/${params.release}/vep_cache/ --stats_file /workdir/filter/${pair_id}.protein -o /workdir/filter/${pair_id}.all_variants.vcf

    filter_vep -i /workdir/filter/${pair_id}.all_variants.vcf --filter "BIOTYPE is protein_coding" | filter_vep --filter "IMPACT is MODERATE" |  filter_vep --filter "IMPACT is HIGH" > /workdir/filter/${pair_id}.protein.vcf
    """
}

process merging {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val pair_id
    val sample_table
    val series

  script:
  """
    #!/usr/local/bin/python

    import pandas as pd
    import numpy as np

    filterfolder="/workdir/filter/"
    sample_sheet="${samplestable}"

    samples = pd.read_excel(sample_sheet, sheet_name='samples', engine='openpyxl')
    samples['SampleID'] = samples['Group'] + '.Rep_' + samples['Replicate'].astype(str)
    sample_dict = samples.set_index('Sample').to_dict()['SampleID']

    # the output file will have the following columns
    # sampleID        chrom   pos     REF     ALT     sums    smaller alternative al. freq. (smaller alternative al.) (%)     genes   ANN     sampleName
    def ANNOTATE(DF, sampleID, sampleName):
      df_format = pd.DataFrame(DF[9].str.split(':').tolist(), columns = ['GT','GQ', 'DP', 'AD', 'VAF', 'PL'])
      DF['sampleID'] = sampleID
      DF['sampleName'] = sampleName
      DF['sums'] = df_format['DP']
      DF['smaller alternative al.'] = df_format['AD'].apply(lambda x: x.split(",")[1] )
      DF['freq. (smaller alternative al.) (%)']= df_format['VAF'].apply(lambda x: x.split(",")[0] ).astype(float) * 100
      DF['Gene'] = DF[7].apply(lambda x: x.split("|")[4])
      DF['SYMBOL'] = DF[7].apply(lambda x: x.split("|")[3])
      DF['Consequence'] = DF[7].apply(lambda x: x.split("|")[1])
      DF['IMPACT'] = DF[7].apply(lambda x: x.split("|")[2])
      DF['Feature_type'] = DF[7].apply(lambda x: x.split("|")[5])
      DF['Feature'] = DF[7].apply(lambda x: x.split("|")[6])
      DF['BIOTYPE'] = DF[7].apply(lambda x: x.split("|")[7])
      DF["ANN"] = DF[7].apply(lambda x: x.split("CSQ=")[1])
      DF=DF[['sampleID', 0,1,3,4,5,"sums","smaller alternative al.", "freq. (smaller alternative al.) (%)", "Gene", "SYMBOL", "Consequence", "IMPACT",\
            "Feature_type", "Feature", "BIOTYPE","ANN", 'sampleName']]
      DF.columns=['sampleID', "chrom","pos","REF","ALT",'QUAL',"sums","smaller alternative al.","freq. (smaller alternative al.) (%)", "Gene",\
                  "SYMBOL", "Consequence", "IMPACT","Feature_type", "Feature", "BIOTYPE","ANN", 'sampleName']
      DF.to_excel(filterfolder + sample_dict[sample['Sample']] + ".results.xlsx", index = False)
      DF.to_csv(filterfolder + sample_dict[sample['Sample']] + ".results.tsv", index = False, sep = '\t')
      return(DF)

    # annotate all samples and merge all samples except wt
    merged_results = pd.DataFrame()
    for index, sample in samples.iterrows():
        print(sample['Sample'])
        if not sample['Sample'] in samples["Background Sample"].tolist():
            tmp = pd.read_csv(filterfolder + sample_dict[sample['Sample']] + ".protein.vcf", sep = '\t', header= None)
            TMP = ANNOTATE(tmp, sample_dict[sample['Sample']], sample['Group'])
            merged_results = merged_results.append(TMP, ignore_index= True)
        else:
            tmp = pd.read_csv(filterfolder + sample_dict[sample['Sample']] + ".protein.vcf", sep = '\t', header= None)
            TMP = ANNOTATE(tmp, sample_dict[sample['Sample']], sample['Group'])
      
    # save master table
    merged_results.to_excel(filterfolder + "${series}.merged_results.xlsx", index = False)
    merged_results.to_csv(filterfolder + "${series}.merged_results.tsv", sep= '\t', index= False)

  """
}

process upload_paths {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
  """
    rm -rf upload.txt

    cd ${params.project_folder}/filter

    for f in \$(ls *results.*) ; do echo "snp_out \$(readlink -f \${f})" >>  upload.txt_ ; done
    
    uniq upload.txt_ upload.txt 
    rm upload.txt_
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
    data = channel.fromFilePairs( "${params.project_folder}/filter/*.SNPs.nowt.vcf", size: -1 )
    vep( data )
    merging( vep.out.collect(), params.samplestable, params.series )
}

workflow upload {
  main:
    upload_paths()
}
