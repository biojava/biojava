package org.biojava.bio.structure.align;

import org.biojava.bio.structure.Atom;

import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.bio.structure.align.model.AFPChain;


public interface StructureAlignment {
	
	/** Run an alignment while specifying the atoms to be aligned. Will used default parameters for the algorithm.
	 * 
	 * @param ca1
	 * @param ca2
	 * @return the afpChain object that contains the alignment.
	 * @throws StructureException
	 */
	public AFPChain align(Atom[] ca1, Atom[] ca2) throws StructureException;
	
	/** run an alignment and also send a bean containing the parameters.
	 * 
	 * @param ca1
	 * @param ca2
	 * @param params
	 * @return the afpChain object that contains the alignment.
	 * @throws StructureException
	 */
	public AFPChain align(Atom[] ca1, Atom[] ca2, Object params) throws StructureException;

	/** Return the paramers for this algorithm. 
	 * 
	 * @return The returned object will be a Java bean.
	 */
	public ConfigStrucAligParams getParameters();
	
	/** Set the default parameters for this algorithm to use
	 * 
	 * @param parameters
	 */
	public void setParameters(ConfigStrucAligParams parameters);
	
	/** Get the name of the Algorithm
	 * 
	 * @return the name of the algorithm
	 */
	public String getAlgorithmName();

	/** Get the Version information for this Algorithm.
	 * 
	 * @return the version of the algorithm
	 */
	public String getVersion();
	
	
	
	/** Returns some documentation on the command line arguments for this algorithm.
	 * 
	 * @return the help string
	 */
	public String printHelp();
	
	
	/** Every alignment algorithm can be called from the command line. The possible arguments are documented in the printHelp method. 
	 */
	//public static void main(String[] argv);
}
