package org.biojava.bio.structure.align;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeSideChainMain;
import org.biojava.bio.structure.align.seq.SmithWaterman3Daligner;


public class StructureAlignmentFactory {

	public static StructureAlignment[] getAllAlgorithms(){
		StructureAlignment[] algorithms = new StructureAlignment[3];
		
		algorithms[0] = new CeMain();
		algorithms[1] = new CeSideChainMain();		
		algorithms[2] = new SmithWaterman3Daligner();
		
		//algorithms[3] = new BioJavaStructureAlignment();
		return algorithms;
		
	}
	
	public static StructureAlignment getAlgorithm(String name) throws StructureException{
		StructureAlignment[] algorithms = getAllAlgorithms();
		for ( StructureAlignment algo : algorithms){
			if (algo.getAlgorithmName().equalsIgnoreCase(name))
				return algo;
		}
		throw new StructureException("Unknown alignment algorithm: " + name);
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
