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

import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;
import org.biojavax.ga.functions.CrossOverFunction;
import org.biojavax.ga.functions.FitnessFunction;
import org.biojavax.ga.functions.MutationFunction;
import org.biojavax.ga.functions.SelectionFunction;

/**
 * The class that runs the cycles of reproduction, evolution and selection,
 * potentially on multiple <code>Population</code>s
 * </p>
 *
 * @author Mark Schreiber
 * @author Susanne Merz
 * @author Andreas Dr&auml;ger
 * @version 1.1
 * @since 1.5
 */

public interface GeneticAlgorithm extends Changeable {

	public static ChangeType	FITNESS_FUNCTION	= new ChangeType(
	                                               "Fitness function changed",
	                                               GeneticAlgorithm.class,
	                                               "FITNESS_FUNCTION");

	/**
	 * The fitness function that will be used to compute the fitness of each
	 * organism.
	 *
	 * @param func
	 *          the <code>FitnessFunction</code> to be used
	 * @throws ChangeVetoException
	 *           if the change is vetoed.
	 */
	public void setFitnessFunction(FitnessFunction func)
	    throws ChangeVetoException;

	/**
	 * Returns the fitness function, i.e. the class that computes the fitness of
	 * each organism in a population.
	 *
	 * @return the fitness function
	 */
	public FitnessFunction getFitnessFunction();

	/**
	 * Sets the <code>Population</code> of <code>Organisms</code> to the
	 * Algorithm.
	 *
	 * @param pop
	 *          the population to add.
	 * @throws ChangeVetoException
	 *           if new populations are not allowed.
	 */
	public void setPopulation(Population pop) throws ChangeVetoException;

	/**
	 * The registered <code>Population</code>
	 *
	 * @return the <code>Population</code> being operated on.
	 */
	public Population getPopulation();

	/**
	 * Changes the <code>SelectionFunction</code> used to select candidates for
	 * the next generation
	 *
	 * @param function
	 *          a <code>SelectionFunction</code>
	 * @throws ChangeVetoException
	 *           if the <code>SelectionFunction</code> is not allowed to be
	 *           changed
	 */
	public void setSelectionFunction(SelectionFunction function)
	    throws ChangeVetoException;

	/**
	 * @return the current <code>SelectionFunction</code>
	 */
	public SelectionFunction getSelectionFunction();

	/**
	 * Changes the <code>CrossOverFunction</code> used to CrossOver Chromosomes
	 *
	 * @param function
	 *          a <code>CrossOverFunction</code>
	 * @throws ChangeVetoException
	 *           if the <code>CrossOverFunction</code> is not allowed to be
	 *           changed
	 */
	public void setCrossOverFunction(CrossOverFunction function)
	    throws ChangeVetoException;

	/**
	 * @return the current CrossOverFunction
	 */
	public CrossOverFunction getCrossOverFunction();

	/**
	 * Sets the current <code>MutationFunction</code>
	 *
	 * @param function
	 *          a <code>MutationFunction</code>
	 * @throws ChangeVetoException
	 *           if the <code>MutationFunction</code> change is Vetoed by a
	 *           listener.
	 */
	public void setMutationFunction(MutationFunction function)
	    throws ChangeVetoException;

	/**
	 * @return the current <code>MutationFunction</code>
	 */
	public MutationFunction getMutationFunction();

	/**
	 * @return the Current generation number
	 */
	public int getGeneration();

	/**
	 * Iterates the Algorithm until the stopping criteria are met. For saftey
	 * implementations should synchronize on this method.
	 *
	 * @param stoppingCriteria
	 *          determines when to stop.
	 * @throws ChangeVetoException
	 *           if the Population being modified is locked
	 * @throws IllegalAlphabetException
	 *           if the MutationFunction chosen attempts to modify a Symbol from
	 *           one of the Chromosomes to a Symbol outside of its Alphabet.
	 * @throws IllegalSymbolException
	 *           if the MutationFunction chosen is using the wrong Alphabet.
	 */
	public void run(GAStoppingCriteria stoppingCriteria)
	    throws ChangeVetoException, IllegalAlphabetException,
	    IllegalSymbolException;

	public static ChangeType	POPULATION	        = new ChangeType(
	                                                  "Population changed",
	                                                  GeneticAlgorithm.class,
	                                                  "POPULATION");

	public static ChangeType	FUNCTION	          = new ChangeType(
	                                                  "Function changed",
	                                                  GeneticAlgorithm.class,
	                                                  "FUNCTION");

	public static ChangeType	CROSS_OVER_FUNCTION	= new ChangeType(
	                                                  "Cross over function changed",
	                                                  GeneticAlgorithm.class,
	                                                  "CROSS_OVER_FUNCTION",
	                                                  FUNCTION);

	public static ChangeType	MUTATION_FUNCTION	  = new ChangeType(
	                                                  "Mutation function changed",
	                                                  GeneticAlgorithm.class,
	                                                  "MUTATION_FUNCTION",
	                                                  FUNCTION);

	public static ChangeType	SELECTION_FUNCTION	= new ChangeType(
	                                                  "Selection function changed",
	                                                  GeneticAlgorithm.class,
	                                                  "SELECTION_FUNCTION",
	                                                  FUNCTION);
}
