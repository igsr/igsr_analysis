/* 
 *  Module containing processes for running the Python codebase from igsr_analysis wrapping the different
 *  vg toolkit commands
 *
 * This module relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 * @author
 * Ernesto Lowy <ernesto.lowy@gmail.com>
 *
 */

 process GIRAFFE {
    /*
    This process run 'vg giraffe'
    */
    tag "$sampleId: $fastq1 $fastq2"
    publishDir "${outdir}", mode: "move", overwrite: true

    maxForks 1
    
    input:
    tuple val(sampleId), val(fastq1), val(fastq2)
    val(gwtname) // optional
    val(graphname) // optional
    val(gbzname) // optional
    val(minimizername)
    val(distname)
    val(outdir)
    val(cpus)

    output:
    path("*.gam")

    script:
    """
    #!/usr/bin/env python
    from GRAPH.VgToolkit import VG
    vg_object = VG()

    outfile = None
    if f"${gbzname}" != "NO_FILE":
        outfile = vg_object.run_giraffe(fastq=f"${fastq1},${fastq2}",
                                        Z=f"${gbzname}",
                                        m=f"${minimizername}",
                                        d=f"${distname}",
                                        t=f"${cpus}",
                                        prefix=f"${sampleId}")
    elif f"${gwtname}" != "NO_FILE":
        outfile = vg_object.run_giraffe(fastq=f"${fastq1},${fastq2}",
                                        H=f"${gwtname}",
                                        g=f"${graphname}",
                                        m=f"${minimizername}",
                                        d=f"${distname}",
                                        t=f"${cpus}",
                                        prefix=f"${sampleId}")
    
    if not outfile:
        raise Exception("No '*.gam' file generated")
    
    print(outfile)
    """
}

process CHUNK {
    /*
    This process run 'vg chunk'
    */
    input:
    val(xgname)
    val(prefix)

    output:
    path("*.pg", emit: chunkFile)

    script:
    """
    #!/usr/bin/env python
    from GRAPH.VgToolkit import VG
    vg_object = VG()
    outfiles = vg_object.run_chunk( n=2,
                                    x=f"${xgname}",
                                    O="pg",
                                    b=f"${prefix}")
    
    if len(outfiles)==0:
        raise Exception("No chunks generated")
    
    for p in outfiles:
        print(p)
    """
}

process AUGMENT {
    /*
    This process run 'vg augment'
    */
    input:
    val(pgname)
    val(gam)

    output:
    path("${pgname.baseName}.aug.pg", emit: aug_pgFile)
    path("${pgname.baseName}.aug.gam", emit: aug_gamFile)

    script:
    """
    #!/usr/bin/env python
    from GRAPH.VgToolkit import VG
    vg_object = VG()
    outfiles = vg_object.run_augment('s',
                                     graph_f=f"${pgname}",
                                     aln_f=f"${gam}",
                                     prefix=f"${pgname.baseName}")

    if len(outfiles)<2:
        raise Exception("No output files")
    
    print(outfiles[0])
    print(outfiles[1])
    """
}

process SNARLS {
    /*
    This process run 'vg snarls'
    */
    input:
    val(graph_f)

    output:
    path("${graph_f.baseName}.snarls", emit: snarlsFile)

    script:
    """
    #!/usr/bin/env python
    from GRAPH.VgToolkit import VG
    vg_object = VG()
    outfile = vg_object.run_snarls(graph_f=f"${graph_f}",
                                   prefix=f"${graph_f.baseName}")
    
    if not outfile:
        raise Exception("No '*.snarls' file generated")
    
    print(outfile)
    """
}

process PACK {
    /*
    This process run 'vg snarls'
    */
    input:
    val(graph_f)
    val(gam)

    output:
    path("${graph_f.baseName}.pack", emit: packFile)

    script:
    """
    #!/usr/bin/env python
    from GRAPH.VgToolkit import VG
    vg_object = VG()
    outfile = vg_object.run_pack(vg_f=f"${graph_f}",
                                 aln_f=f"${gam}",
                                 prefix=f"${graph_f.baseName}",
                                 Q=5)
    
    if not outfile:
        raise Exception("No '*.pack' file generated")
    
    print(outfile)
    """
}

process CALL {
    /*
    This process run 'vg call'
    */
    input:
    val(graph_f)
    val(snarl_f)
    val(pack_f)

    output:
    path("${graph_f.baseName}.vcf", emit: vcfFile)

    script:
    """
    #!/usr/bin/env python
    from GRAPH.VgToolkit import VG
    vg_object = VG()
    outfile = vg_object.run_call(graph_f=f"${graph_f}",
                                 prefix=f"${graph_f.baseName}",
                                 k=f"${pack_f}",
                                 r=f"${snarl_f}")
    
    if not outfile:
        raise Exception("No '*.vcf' file generated")
    
    print(outfile)
    """
}

