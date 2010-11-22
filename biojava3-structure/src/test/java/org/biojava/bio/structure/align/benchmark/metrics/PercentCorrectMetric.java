package org.biojava.bio.structure.align.benchmark.metrics;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.benchmark.MultipleAlignment;
import org.biojava.bio.structure.align.model.AFPChain;

/**
 * Calculates the percentage of the reference alignment which was scored correctly:<br/>
 * (# residue pairs in the reference alignment which appear in the test alignment) /
 * (# pairs in the reference alignment)
 * <p>
 * Note that defining additional pairs is not detrimental, so this is appropriate
 * for limited manual alignments, eg of active sites.
 * @author Spencer Bliven
 *
 */
public class PercentCorrectMetric extends Metric {

	@Override
	public double calculate(MultipleAlignment reference, AFPChain align, Atom[] ca1, Atom[] ca2, Map<String, Object> metaData) {
		
		List<Atom[]> structures = new ArrayList<Atom[]>(2);
		structures.add(ca1);
		structures.add(ca2);
		int[][] refAln;
		try {
			refAln = reference.getAlignmentMatrix(structures);
		} catch (StructureException e) {
			e.printStackTrace();
			return Double.NaN;
		}

		int[][][] optAln = align.getOptAln();
		int[] blockLens = align.getOptLen();

		// Hash table of ca1->ca2 values for the alignment
		HashMap<Integer,Integer> optAlnMap = new HashMap<Integer,Integer>();
		for(int block=0;block<align.getBlockNum(); block++) {
			for(int i=0;i<blockLens[block];i++) {
				optAlnMap.put(optAln[block][0][i], optAln[block][1][i]);
			}
		}
		
		// Calculate number correct
		int correct = 0;
		for(int i=0;i<reference.size();i++) {
			if(optAlnMap.containsKey(refAln[0][i])) {
				if(optAlnMap.get(refAln[0][i]).equals(refAln[1][i])) {
					correct++;
				}
			}
		}
		
		return 100.* correct/reference.size();
	}

	@Override
	public String getName() {
		return "%_correct";
	}
	
	@Override
	public String format(double result) {
		return String.format("%02.2f",result);
	}

}
