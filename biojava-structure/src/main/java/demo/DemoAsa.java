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
package demo;



import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.asa.AsaCalculator;
import org.biojava.nbio.structure.asa.GroupAsa;

import java.io.IOException;

public class DemoAsa {

	private static final boolean hetAtoms = false;

	public static void main(String[] args) throws IOException, StructureException {

		String pdbCode = args[0];
		int numThreads = Integer.parseInt(args[1]);

		demoAsa(pdbCode, numThreads);
	}

	private static void demoAsa(String pdbCode, int numThreads) throws IOException, StructureException {

		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);

		StructureIO.setAtomCache(cache);

		Structure structure = StructureIO.getStructure(pdbCode);

		long start = System.currentTimeMillis();

		AsaCalculator asaCalc = new AsaCalculator(structure,
				AsaCalculator.DEFAULT_PROBE_SIZE,
				1000, numThreads, hetAtoms);

		GroupAsa[] groupAsas = asaCalc.getGroupAsas();

		long end = System.currentTimeMillis();


		double tot = 0;



		for (GroupAsa groupAsa: groupAsas) {
			System.out.printf("%1s\t%5s\t%3s\t%6.2f\n",
					groupAsa.getGroup().getChainId(),
					groupAsa.getGroup().getResidueNumber(),
					groupAsa.getGroup().getPDBName(),
					groupAsa.getAsaU());
			tot+=groupAsa.getAsaU();
		}


		System.out.printf("Total area: %9.2f\n",tot);
		System.out.printf("Time: %4.1fs\n",((end-start)/1000.0));


		System.out.println("Testing scaling: ");
		double[] runTimes = new double[numThreads];
		for (int nThreads=1;nThreads<=numThreads;nThreads++) {
			start = System.currentTimeMillis();

			asaCalc = new AsaCalculator(structure,
					AsaCalculator.DEFAULT_PROBE_SIZE,
					1000, numThreads, hetAtoms);

			// only calculating all atom ASAs without keeping the returned value
			asaCalc.calculateAsas();

			end = System.currentTimeMillis();
			runTimes[nThreads-1] = (end-start)/1000.0;

		}
		for (int nThreads=1;nThreads<=numThreads;nThreads++) {
			System.out.printf(nThreads+" threads, time: %4.1fs -- x%2.1f\n",runTimes[nThreads-1],runTimes[0]/runTimes[nThreads-1]);
		}


	}


}
