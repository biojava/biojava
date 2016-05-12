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

import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;

/**
 * Interface for the Multiple Structure Alignment Algorithms. The Pairwise Alignment Algorithms can also
 * implement this class to be able to return {@link MultipleAlignment} Objects.
 *
 * @author Aleix Lafita
 *
 */
public interface MultipleStructureAligner{

	/**
	 * Run an alignment while specifying the atoms to be aligned.
	 * The default parameters for the algorithm are used.
	 *
	 * @param atomArrays List of Atoms of all the structures
	 * @return MultipleAlignment object that contains the alignment.
	 * @throws StructureException
	 * @see #align(List,Object)
	 */
	public MultipleAlignment align(List<Atom[]> atomArrays) throws StructureException;

	/**
	 * Run an alignment and also send a bean containing the parameters.
	 *
	 * @param atomArrays List of Atoms of all the structures
	 * @return MultipleAlignment object that contains the alignment.
	 * @throws StructureException
	 * @see #align(List)
	 */
	public MultipleAlignment align(List<Atom[]> atomArrays, Object params) throws StructureException;

	/**
	 * Return the parameters of this algorithm instance.
	 *
	 * @return The returned Object will be a Java bean.
	 */
	public ConfigStrucAligParams getParameters();

	/**
	 * Set the parameters for this algorithm to use.
	 *
	 * @param parameters
	 */
	public void setParameters(ConfigStrucAligParams parameters);

	/**
	 * Get the name of this Algorithm.
	 *
	 * @return String name of the algorithm
	 */
	public String getAlgorithmName();

	/**
	 * Get the Version information for this Algorithm.
	 *
	 * @return String version of the algorithm
	 */
	public String getVersion();
}
