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
    input:
    val(fastq)
    path(gbwtFile)
    path(graphname)
    path(minf)
    path(distFile)
    val(cpus)
    val(prefix)

    script:
    """
    #!/usr/bin/env python
    from GRAPH.VgToolkit import VG
    vg_object = VG()
    outfile = vg_object.run_giraffe(fastq=f"${fastq}",
                                    H=f"${gbwtFile}",
                                    g=f"${graphname}",
                                    m=f"${minf}",
                                    d=f"${distFile}",
                                    t=f"${cpus}",
                                    prefix=f"${prefix}")
    """
}


