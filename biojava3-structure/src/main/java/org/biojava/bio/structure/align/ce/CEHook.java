/**
 * 
 */
package org.biojava.bio.structure.align.ce;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.model.AFPChain;

/**
 * The CE algorithm allows generic alignment of structures. This can be used as
 * a component of modified algorithms, such as CE-CP, by modifying the input 
 * structures and then fixing up the output. While this could be done using a
 * custom @{link StructureAlignment} algorithm, it is sometimes useful to do
 * this from within CeMain based on parameter settings. This class provides
 * hooks to modify structures before and after alignment.
 * 
 * @author Spencer Bliven
 *
 */
public interface CEHook {

	/**
	 * Perform preprocessing of Atom arrays before running CE.
	 * @param ca1 Original CA atoms of first structure
	 * @param ca2 Original CA atoms of second structure
	 * @param param A parameter object, typically a CeParameters instance
	 * @return A length-2 array giving modified {ca1, ca2}
	 */
	public Atom[][] preprocess(Atom[] ca1, Atom[] ca2, Object param);
	/**
	 * Perform postprocessing of Atom arrays before running CE.
	 * @param ca1m Modified CA atoms of first structure, after preprocessing
	 * @param ca2m Modified CA atoms of second structure, after preprocessing
	 * @param param A parameter object, typically a CeParameters instance
	 * @param afpChain AFPChain giving the alignment between ca1m and ca2m
	 * @return A modified AFPChain giving the alignment between the original ca1 and ca2
	 */
	public AFPChain postprocess(Atom[] ca1m, Atom[] ca2m, Object param, AFPChain afpChain);
	
}
