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
import org.biojava.utils.ChangeVetoException;
import org.biojavax.ga.GeneticAlgorithm;
import org.biojavax.ga.Organism;
import org.biojavax.ga.Population;
import org.biojavax.ga.exception.IllegalOrganismException;
import org.biojavax.ga.util.WeightedSet;

/**
 * <p>
 * A Selection function that determines the proportion of individuals in a new
 * population proportionally to their fitness. The population size is not
 * allowed to grow. Individuals are randomly selected for replication, those
 * with greater fitness tend to replicate more often.
 *
 * @author Mark Schreiber
 * @author Susanne Merz
 * @author Andreas Dr&auml;ger
 * @version 1.1
 */

public class ProportionalSelection implements SelectionFunction {

	public ProportionalSelection() {

	}

	public Population select(Population pop, GeneticAlgorithm genAlg)
	    throws ChangeVetoException {
		WeightedSet set = new WeightedSet();
		int size = pop.size();

		for (Iterator i = pop.organisms(); i.hasNext();) {
			Object item = i.next();
			double fit[] = ((Organism) item).getFitness();
			// TODO: maybe we have to consider every fitness value.
			set.setWeight(item, fit[0]);
		}

		pop.removeAllOrganisms();

		for (int i = 0; i < size; i++) {
			try {
				Organism o = (Organism) set.sample();
				String name = o.getName().split(":")[0];

				// to begin with the name may not have a ":" in it
				if (name.equals("")) {
					name = o.getName();
				}

				// System.out.println("name: "+name);
				o = o.replicate(name + ":" + i);
				pop.addOrganism(o);
			} catch (IllegalOrganismException ex) {
				throw new BioError("A previously legal organism is now illegal??", ex);
			}
		}

		return pop;
	}
}
