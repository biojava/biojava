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
 * A metric based off the concept of Alignment consistency, as presented in<br/>
 * <i>Mayr et al. Comparative analysis of protein structure alignments. BMC Struct Biol (2007) vol. 7 pp. 50</i>.
 * <p>
 * The formula is A(s) = I(s)/L_ref, where A(s) is the alignment consistency, 
 * I(s) is the number of pairs in the reference alignment for which there exists
 * some pair in the other alignment which matches within s residues.
 * <p>
 * The metric is slightly modified from that in the paper; rather than using the
 * longer alignment length as the denominator we use the reference alignment
 * length. Thus A(0) is equivalent to {@link PercentCorrectMetric}/100.
 * @author Spencer Bliven
 *
 */
public class ConsistencyMetric extends Metric {
	int shift; // 's' in the paper

	public ConsistencyMetric() { this(0); }
	public ConsistencyMetric(int s) {
		this.shift = s;
	}

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
		HashMap<Integer,Integer> optAlnMap1 = new HashMap<Integer,Integer>();
		// Hash table of ca2->ca1 values for the alignment
		HashMap<Integer,Integer> optAlnMap2 = new HashMap<Integer,Integer>();
		for(int block=0;block<align.getBlockNum(); block++) {
			for(int i=0;i<blockLens[block];i++) {
				optAlnMap1.put(optAln[block][0][i], optAln[block][1][i]);
				optAlnMap2.put(optAln[block][1][i], optAln[block][0][i]);
			}
		}

		// Calculate number correct
		int correct = 0;
		for(int i=0;i<reference.size();i++) {
			int ca1ref = refAln[0][i];
			int ca2ref = refAln[1][i];
			if(optAlnMap1.containsKey(ca1ref) ) {
				if( Math.abs(optAlnMap1.get(ca1ref).intValue() - ca2ref) <= shift ) {
					correct++;
				}
			} else if(optAlnMap2.containsKey(ca2ref) ) {
				if( Math.abs(optAlnMap2.get(ca2ref).intValue() - ca1ref) <= shift ) {
					correct++;
				}
			}
		}

		return (double)correct/reference.size();
	}

	@Override
	public String getName() {
		return "A_"+shift;
	}

	@Override
	public String format(double result) {
		return String.format("%.4f",result);
	}
	
	
	/**
	 * @return the shift
	 */
	public int getShift() {
		return shift;
	}
	/**
	 * @param s the shift to set
	 */
	public void setShift(int s) {
		this.shift = s;
	}

}
