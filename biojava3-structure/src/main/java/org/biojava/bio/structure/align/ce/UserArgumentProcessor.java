package org.biojava.bio.structure.align.ce;

public interface UserArgumentProcessor {

	
	/** Process user arguments that have been provided from the command line
	 * 
	 * @param argv
	 */
	public void process(String[] argv);
	
	/**
	 * Print help about the arguments
	 * @return
	 */
	public String printHelp();
}
