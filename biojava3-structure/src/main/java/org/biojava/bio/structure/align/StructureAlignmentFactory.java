package org.biojava.bio.structure.align;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.fatcat.FatCatFlexible;
import org.biojava.bio.structure.align.fatcat.FatCatRigid;
import org.biojava.bio.structure.align.seq.SmithWaterman3Daligner;


public class StructureAlignmentFactory {

	private static List<StructureAlignment> algorithms = new ArrayList<StructureAlignment>();
	
	static {
		algorithms.add( new CeMain() );
		algorithms.add( new CeCPMain() );
		//algorithms.add(new CeSideChainMain());

		StructureAlignment fatcatRigid    = new FatCatRigid();
		StructureAlignment fatcatFlexible = new FatCatFlexible();

		if ( fatcatRigid != null) {
			algorithms.add( fatcatRigid) ;

		}
		if ( fatcatFlexible != null){
			algorithms.add( fatcatFlexible );

		}
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

	public static StructureAlignment getAlgorithm(String name) throws StructureException{
		StructureAlignment[] algorithms = getAllAlgorithms();
		for ( StructureAlignment algo : algorithms){
			if (algo.getAlgorithmName().equalsIgnoreCase(name))
				return algo;
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

		return (String[])names.toArray(new String[names.size()]);
	}

}
