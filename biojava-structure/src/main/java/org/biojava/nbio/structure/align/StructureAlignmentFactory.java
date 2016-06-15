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

import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.ce.CeCPMain;
import org.biojava.nbio.structure.align.ce.CeMain;
import org.biojava.nbio.structure.align.fatcat.FatCatFlexible;
import org.biojava.nbio.structure.align.fatcat.FatCatRigid;
import org.biojava.nbio.structure.align.seq.SmithWaterman3Daligner;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;


public class StructureAlignmentFactory {

	private final static Logger logger = LoggerFactory.getLogger(StructureAlignmentFactory.class);

	private static List<StructureAlignment> algorithms = new ArrayList<StructureAlignment>();

	static {
		algorithms.add( new CeMain() );
		algorithms.add( new CeCPMain() );
		//algorithms.add( new OptimalCECPMain() );
		//algorithms.add(new CeSideChainMain());

		StructureAlignment fatcatRigid    = new FatCatRigid();
		StructureAlignment fatcatFlexible = new FatCatFlexible();

		algorithms.add( fatcatRigid) ;

		algorithms.add( fatcatFlexible );

		algorithms.add( new SmithWaterman3Daligner()) ;
		//algorithms.add( new BioJavaStructureAlignment());
	}

	/**
	 * Adds a new StructureAlignment algorithm to the list.
	 *
	 * Only one instance is stored for each algorithmName, so it is possible
	 * that a different instance may be returned by getAlgorithm(alg.getAlgorithmName())
	 *
	 * @param alg the alignment algorithm
	 */
	public static void addAlgorithm(StructureAlignment alg) {
		//ensure uniqueness
		try {
			getAlgorithm(alg.getAlgorithmName());
			// algorithm was found. Do nothing.
		} catch(StructureException e) {
			// no algorithm found, so it's new
			algorithms.add(alg);
		}
	}

	/**
	 * Removes the specified algorithm from the list of options
	 * @param name the name of the algorithm to remove
	 * @return true if the specified algorithm was found and removed
	 */
	public static boolean removeAlgorithm(String name) {
		ListIterator<StructureAlignment> algIt = algorithms.listIterator();
		while(algIt.hasNext()) {
			StructureAlignment alg = algIt.next();
			if(alg.getAlgorithmName().equalsIgnoreCase(name)) {
				algIt.remove();
				return true;
			}
		}
		return false;
	}

	/**
	 * Removes all algorithms from the list
	 */
	public static void clearAlgorithms() {
		algorithms.clear();
	}

	public static StructureAlignment getAlgorithm(String name) throws StructureException{
		for ( StructureAlignment algo : algorithms){
			if (algo.getAlgorithmName().equalsIgnoreCase(name)) {
				//return algo;
				// CeCalculator is not thread safe,
				// avoid issues with this in multi-threaded environments bu
				// creating a new StructureAlignment every time this is called
				try {
					@SuppressWarnings("unchecked")
					Class<StructureAlignment> c = (Class<StructureAlignment>) Class.forName(algo.getClass().getName());
					return c.newInstance();
				} catch (ClassNotFoundException e){
					logger.error("Exception: ", e);
					return null;
				} catch (IllegalAccessException e){
					logger.error("Exception: ", e);
					return null;
				} catch (InstantiationException e){
					logger.error("Exception: ", e);
					return null;
				}


			}
		}

		throw new StructureException("Unknown alignment algorithm: " + name);
	}

	public static StructureAlignment[] getAllAlgorithms(){
		return algorithms.toArray(new StructureAlignment[algorithms.size()]);
	}

	public static String[] getAllAlgorithmNames(){
		StructureAlignment[] algos = getAllAlgorithms();
		List<String> names = new ArrayList<String>();

		for (StructureAlignment alg : algos){
			names.add(alg.getAlgorithmName());
		}

		return names.toArray(new String[names.size()]);
	}

}
