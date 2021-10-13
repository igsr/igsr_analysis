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

    Parameters
    ----------
    fastq : String with a fastq or comma-separated fastq paths

    Returns
    -------
    TODO
    */
    tag "$sampleId: $fastq1 $fastq2"

    input:
    tuple val(sampleId), val(fastq1), val(fastq2)
    val(gwtname)
    val(graphname)
    val(minimizername)
    val(distname)
    val(cpus)

    script:
    """
    #!/usr/bin/env python
    from GRAPH.VgToolkit import VG
    vg_object = VG()
    outfile = vg_object.run_giraffe(fastq=f"${fastq1},${fastq2}",
                                    H=f"${gwtname}",
                                    g=f"${graphname}",
                                    m=f"${minimizername}",
                                    d=f"${distname}",
                                    t=f"${cpus}",
                                    prefix=f"${sampleId}")
    """
}

process VG_CALL {
    
}


