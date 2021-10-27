process SAVE_FILE {
    /*
    This process is of general utility and used to save a file
	in a directory
	
	Parameters
	----------
	file : path to file to save
	dirname : name of directory used to save 'file'
	prefix : output name for saved file
	mode : mode used to save the file: ['move','copy']
	*/
	publishDir "${dirname}", mode: "${mode}", overwrite: true

    input:
		path(afile)
	    val(dirname)
	    val(prefix)
		val(mode)

    output:
        path "${prefix}"

    """
    mv ${afile} ${prefix}
    """
    }
    