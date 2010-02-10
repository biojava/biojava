/*
 * BioJava development code This code may be freely distributed and modified
 * under the terms of the GNU Lesser General Public Licence. This should be
 * distributed with the code. If you do not have a copy, see:
 * http://www.gnu.org/copyleft/lesser.html Copyright for this code is held
 * jointly by the individual authors. These should be listed in @author doc
 * comments. For more information on the BioJava project and its aims, or to
 * join the biojava-l mailing list, visit the home page at:
 * http://www.biojava.org/
 */

package org.biojavax.ga;

import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;

/**
 * A GA 'organism' contains one or more Chromosomes
 *
 * @author Mark Schreiber
 * @author Susanne Merz
 * @author Andreas Dr&auml;ger
 * @version 1.0
 * @since 1.5
 */

public interface Organism extends Changeable {

	/**
	 * This method allows to set the fitness of this organism to the specified
	 * value. Generally this will be an array, which in the most cases contains
	 * just a single entry. In cases where we want to have multi-objective
	 * optimization we may want to make use of a more general fitness array with
	 * mutliple entries.
	 *
	 * @param fitness
	 */
	public void setFitness(double[] fitness);

	/**
	 * Returns the current fitness of this organism. This is an array. Note that
	 * in the most cases this array may only contain one single value, but for
	 * multi-objective optimization it is necessary to store multiple fitness
	 * values.
	 *
	 * @return the fitness of the organism
	 */
	public double[] getFitness();

	/**
	 * Gets the organisms 'chromosome' sequences
	 *
	 * @return a <code>SymbolList[]</code>
	 */
	public SymbolList[] getChromosomes();

	/**
	 * Sets the organisms 'chromosome' sequences.
	 *
	 * @param chromosomes
	 *          a <code>SymbolList[]</code>
	 * @throws ChangeVetoException
	 *           if the Chromosome collection of the Organism is unchangable
	 */
	public void setChromosomes(SymbolList[] chromosomes)
	    throws ChangeVetoException;

	/**
	 * Gets the organisms name
	 *
	 * @return the name String
	 */
	public String getName();

	/**
	 * Sets the organisms name
	 *
	 * @param name
	 *          the name of the organism.
	 * @throws ChangeVetoException
	 *           if the name may not be changed.
	 */
	public void setName(String name) throws ChangeVetoException;

	/**
	 * Creates a replica of this <code>Organism</code> with a new name.
	 *
	 * @param name
	 *          the new name for the sequence.
	 * @return the replicated organism.
	 */
	public Organism replicate(String name);

	/**
	 * Is the organism Haploid?
	 *
	 * @return true if it is.
	 */
	public boolean isHaploid();

	public static final ChangeType	CHROMOSOMES	= new ChangeType(
	                                                "Chromosomes changed",
	                                                "ga.Organism", "CHROMOSOMES");

	public static final ChangeType	NAME	      = new ChangeType("Name changed",
	                                                "ga.Organism", "NAME");

}
