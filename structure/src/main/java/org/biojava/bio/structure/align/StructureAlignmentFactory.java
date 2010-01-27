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

		// check if the fatcat source code is in the class path
		// we still need to seek permission to get jFatCat released under the LGPL
		if  ( name.equals("jFatCat_rigid")){
			return getFatCatRigid();
		}
		
		if ( name.equals("jFatCat_flexible")){
			return getFatCatFlexible();
		}
		throw new StructureException("Unknown alignment algorithm: " + name);
	}

	private static StructureAlignment getFatCatRigid() throws StructureException{

		String className = "org.rcsb.fatcat.FatCatRigid";
		try {
			Class c = Class.forName(className);
			Object o = c.newInstance();
			return (StructureAlignment)o;
		} catch (Exception e){
			throw new StructureException("Could not find FatCatRigid in the classpath.");
		}
	}
	
	private static StructureAlignment getFatCatFlexible() throws StructureException{

		String className = "org.rcsb.fatcat.FatCatFlexible";
		try {
			Class c = Class.forName(className);
			Object o = c.newInstance();
			return (StructureAlignment)o;
		} catch (Exception e){
			throw new StructureException("Could not find FatCatFlexible in the classpath.");
		}
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
