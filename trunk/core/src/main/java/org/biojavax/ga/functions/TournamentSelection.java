package org.biojavax.ga.functions;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.Random;

import org.biojava.bio.BioError;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.ga.GeneticAlgorithm;
import org.biojavax.ga.Organism;
import org.biojavax.ga.Population;
import org.biojavax.ga.exception.IllegalOrganismException;
import org.biojavax.ga.functions.SelectionFunction;
import org.biojavax.ga.impl.SimplePopulation;

/**
 * Tournament Selection chooses the best organisms from n random subsets of a
 * given population. Currently it assumes a maximization problem. Perhaps this
 * could be selected depending on the Genetic Algorithm utilized.
 *
 * @author Susanne Merz
 */

public class TournamentSelection implements SelectionFunction {

	private int	selpressure;

	/**
	 * Default constructor: sets the selection pressure to the value of 10.
	 */
	public TournamentSelection() {
		this.setSelectionPressure(10);
	};

	/**
	 * sets the parameter controlling selection pressure
	 *
	 * @param numberOfIndividuals
	 *          the number of Individuals the best is selected from, ranges from 1
	 *          (random selection) to the size of the population (elitism)
	 */
	public void setSelectionPressure(int numberOfIndividuals) {
		selpressure = numberOfIndividuals;
	}

	/**
	 * Standard call to select organisms, will select a number of Organisms
	 * corresponding to 75 % of the population.
	 *
	 * @see org.biojavax.ga.functions.SelectionFunction#select(org.biojavax.ga.Population, org.biojavax.ga.GeneticAlgorithm)
	 */
	public Population select(Population pop, GeneticAlgorithm genAlg)
	    throws ChangeVetoException {
		return selectNIndividuals(pop, genAlg, pop.size());

	}

	/**
	 * This method selects n Organism from the population it is given, using the
	 * tournament selection method
	 *
	 * @param pop
	 *          the population to select from
	 * @param ga
	 *          the <code>GeneticAlgorithm</code> this selection belongs to
	 * @param n
	 *          number of individuals to be selected.
	 * @return nextgen a <code>Population</code> containing the selected
	 *         Organisms
	 */
	public Population selectNIndividuals(Population pop, GeneticAlgorithm ga,
	    int n) {

		SimplePopulation nextgen = new SimplePopulation();
		if (selpressure <= 0) {
			System.out.println("Sorry can't select with a selection pressure <= 0");
			return null;
		}
		int q = selpressure;

		Random r = new Random();
		// choose n individuals
		for (int i = 0; i < n; i++) {
			// need to check that the size of the set we want
			// to choose from does not exceed the populations size
			if (q > pop.size()) {
				q = pop.size();
			}
			// Maximization Problem
			double bestvalue = Double.MIN_VALUE;
			Organism best = null;
			// create a List of Organisms, simulating a subpopulation of size q
			LinkedList<Organism> subpopulation = new LinkedList<Organism>();
			while (subpopulation.size() < q) {
				Iterator<Organism> it = pop.organisms();
				int pos = r.nextInt(pop.size());
				Organism current = it.next();
				while (pos > 0) {
					current = it.next();
					pos--;
				}
				// if the subpopulation already contains the choosen
				// element, we pick the next that isn't contained
				while (subpopulation.contains(current)) {
					// run in circles ;)
					if (!it.hasNext()) {
						it = pop.organisms();
					}
					current = it.next();
				}
				// while we already have a handle on the organism
				// we want save which one is the best so far in the
				// subpopulation
				double value[] = current.getFitness();
				if ((value[0] >= bestvalue)) {
					best = current;
					bestvalue = value[0];
				}

				subpopulation.add(current);

			}

			// finished creating subpopulation, add best organism
			// from the subpopulation to the next generation

			try {
				if (best != null) {
					String namesprefix = best.getName().split(":")[0];
					Organism o = best.replicate(namesprefix + ":" + i);
					nextgen.addOrganism(o);
					// pop.removeOrganism(best);
				}
			} catch (IllegalOrganismException e) {
				e.printStackTrace();
			}
		}
		pop.removeAllOrganisms();
		for (Iterator it = nextgen.organisms(); it.hasNext();) {
			try {
				pop.addOrganism((Organism) it.next());
			} catch (ChangeVetoException e) {

				e.printStackTrace();
			} catch (IllegalOrganismException ex) {
				throw new BioError("A previously legal organism is now illegal??", ex);
			}

		}
		pop = nextgen;
		return nextgen;
	}

}
