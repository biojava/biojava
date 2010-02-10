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

package org.biojavax.ga.functions;

import org.biojavax.ga.GeneticAlgorithm;
import org.biojavax.ga.Organism;
import org.biojavax.ga.Population;

/**
 * Calculates the fitness of an <code>Organism</code> in a
 * <code>Population</code> of <code>Organisms</code>
 *
 * @author Mark Schreiber
 * @author Susanne Merz
 * @author Andreas Dr&auml;ger
 * @version 1.1
 * @since 1.5
 */

public interface FitnessFunction {

	/**
	 * Calculates the fitness of <code>org</code>. This can be done
	 * independently of the Population pop (by ignoring the argument in your
	 * implementation) or dependent on the other members of the
	 * <code>Population pop</code>. Every implementation of this function
	 * should assign the fitness value computed in this function to the given
	 * organism. This is important so that the organism knows its current fitness.
	 * Note that for simple problems this fitness array will contain only one
	 * single value. However, to enable multi-objective fitness functions this has
	 * to be an array.
	 *
	 * @param org
	 *          The <code>Organism</code> to score
	 * @param pop
	 *          The <code>Population</code> to consider
	 * @param genAlg
	 *          the parent<code>GeneticAlgorithm</code>
	 * @return the fitness score.
	 */
	public double[] fitness(Organism org, Population pop, GeneticAlgorithm genAlg);
}
