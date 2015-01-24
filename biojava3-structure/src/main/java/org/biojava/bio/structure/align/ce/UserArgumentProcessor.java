package org.biojava.bio.structure.align.ce;

import org.biojava.bio.structure.align.util.ConfigurationException;

public interface UserArgumentProcessor {

	
	/** Process user arguments that have been provided from the command line
	 * 
	 * @param argv
	 */
	public void process(String[] argv) throws ConfigurationException;
	
	/**
	 * Print help about the arguments
	 * @return
	 */
	public String printHelp();
}
