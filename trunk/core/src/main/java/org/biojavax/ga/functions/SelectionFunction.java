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

import java.util.Iterator;

import org.biojava.bio.BioError;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.ga.GeneticAlgorithm;
import org.biojavax.ga.Organism;
import org.biojavax.ga.Population;

/**
 * Selects Organisms for Replication and returns the offspring.
 *
 * @author Mark Schreiber
 * @author Susanne Merz
 * @author Andreas Dr&auml;ger
 * @version 1.1
 * @since 1.5
 */

public interface SelectionFunction {
	/**
	 * Selects all members of a population for replication
	 */
	public static SelectionFunction	DEFAULT	= new SelectAll();

	/**
	 * Selects a <code>Population</code> of <code>Organisms</code> for
	 * replication based on their fitness.
	 *
	 * @param pop
	 *          the <code>Population</code> to select from.
	 * @param genAlg
	 *          the parent <code>GeneticAlgorithm</code>.
	 * @return the <code>Organism</code>s selected
	 * @throws ChangeVetoException
	 *           if the function attempts to change the population and it is
	 *           vetoed.
	 */
	public Population select(Population pop, GeneticAlgorithm genAlg)
	    throws ChangeVetoException;

	/*-----------INNER CLASSES---------------------------------*/

	/**
	 * <p>
	 * Selects individuals who's fitness exceeds a threshold value.
	 * </p>
	 *
	 * @author Mark Schreiber
	 * @version 1.0
	 */
	public final class Threshold implements SelectionFunction {
		private double	cutoff;

		public Threshold(double cutoff) {
			this.cutoff = cutoff;
		}

		public double getCutoff() {
			return cutoff;
		}

		/**
		 * Selects individuals whose fitness (as determined by the
		 * <code>FitnessFunction</code>) is more than the cutoff. Removes those
		 * that aren't.
		 *
		 * @param pop
		 *          the <code>Population</code> to select from.
		 * @param genAlg
		 *          the parent <code>GeneticAlgorithm</code>
		 * @return the <code>Population</code> of selected individuals.
		 */
		public Population select(Population pop, GeneticAlgorithm genAlg) {
			for (Iterator i = pop.getOrganisms().iterator(); i.hasNext();) {
				Organism o = (Organism) i.next();
				try {
					double fitness[] = o.getFitness();
					boolean remove = false;
					for (int j = 0; (j < fitness.length) && !remove; j++)
						if (fitness[j] < cutoff) remove = true;
					if (remove) {
						pop.removeOrganism(o);
						// System.out.println("removing organism "+o.getName());
					}
				} catch (ChangeVetoException ex) {
					throw new BioError(
					    "population has been locked, cannot select individuals", ex);
				}
			}
			return pop;
		}
	}

	public final class SelectAll implements SelectionFunction {
		public Population select(Population pop, GeneticAlgorithm genAlg) {
			return pop;
		}

		/**
		 * @throws UnsupportedOperationException
		 *           as there is no <code>FitnessFunction</code> for this class
		 * @return you won't get this far, trust me!
		 */
		public FitnessFunction getFitnessFunction() {
			throw new UnsupportedOperationException(
			    "No FitnessFunction defined for SelectAll SelectionFunction");
		}

		/**
		 * @param func
		 *          you could try this but it will throw a
		 *          <code>ChangeVetoException</code>
		 * @throws ChangeVetoException
		 *           you Cannot set the <code>FitnessFunction</code> of a
		 *           <code>SelectAll</code> <code>SelectionFunction</code>"
		 */
		public void setFitnessFunction(FitnessFunction func)
		    throws ChangeVetoException {
			throw new ChangeVetoException(
			    "Cannot set the FitnessFunction of SelectAll SelectionFunction");
		}
	}
}
