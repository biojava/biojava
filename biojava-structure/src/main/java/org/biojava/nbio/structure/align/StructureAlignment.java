/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.structure.align;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.nbio.structure.align.model.AFPChain;


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
}
