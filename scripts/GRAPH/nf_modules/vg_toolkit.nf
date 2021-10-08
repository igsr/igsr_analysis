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
    fastq : FASTQ file to be aligned

    Returns
    -------
    TODO
    */
    input:
    path(fastq)
    path(gbwtf)
    path(gbwtgraphf)
    path(minf)
    path(distf)

    script:
    """
    #!/usr/bin/env python
    from GRAPH.VgToolkit import VG
    print("${fastq}")
    vg_object = VG()
    outfile = vg_object.run_giraffe(gbz_f=f"{datadir}/outdir/test.autoindex.giraffe.gbz",
                                    min=f"{datadir}/outdir/test.autoindex.min",
                                    dist=f"{datadir}/outdir/test.autoindex.dist",
                                    fastq=f"{datadir}/VG/input.fq",
                                    prefix=f"{datadir}/outdir/test")
    """
}


