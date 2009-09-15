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

import junit.framework.TestCase;

import org.biojava.utils.ChangeVetoException;
import org.biojavax.ga.GeneticAlgorithm;
import org.biojavax.ga.Organism;
import org.biojavax.ga.Population;
import org.biojavax.ga.impl.SimpleGeneticAlgorithm;
import org.biojavax.ga.impl.SimpleOrganism;
import org.biojavax.ga.impl.SimplePopulation;

/**
 * @author Mark Schreiber
 */
public class ProportionalSelectionTest extends TestCase {
	private ProportionalSelection	ps;

	private Population	          pop;

	private GeneticAlgorithm	    ga;

	public ProportionalSelectionTest(String s) {
		super(s);
	}

	protected void setUp() throws Exception {
		super.setUp();

		pop = new SimplePopulation("test");

		for (int i = 0; i < 5; i++) {
			pop.addOrganism(new SimpleOrganism("a" + i));
			pop.addOrganism(new SimpleOrganism("b" + i));
		}

		ps = new ProportionalSelection();
		ga = new SimpleGeneticAlgorithm(pop, MutationFunction.NO_MUTATION,
		    CrossOverFunction.NO_CROSS, ps);
		ga.setFitnessFunction(new FitnessFunction() {
			public double[] fitness(Organism org, Population pop,
			    GeneticAlgorithm genAlg) {
				return (org.getName().startsWith("a")) ? new double[] {1.0}
				    : new double[] {0.0};
			}
		});
	}

	protected void tearDown() throws Exception {
		ps = null;
		pop = null;
		ga = null;
		super.tearDown();
	}

	public void testSelect() {
		try {
			ps.select(pop, ga);
		} catch (ChangeVetoException ex) {
			fail(ex.getMessage());
		}

		int count = 0;
		for (Iterator i = pop.organisms(); i.hasNext();) {
			Organism o = (Organism) i.next();
			if (o.getName().startsWith("a")) count++;
		}
		assertTrue(count == 10);
	}
}
