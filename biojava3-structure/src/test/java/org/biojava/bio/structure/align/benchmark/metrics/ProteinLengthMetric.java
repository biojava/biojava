package org.biojava.bio.structure.align.benchmark.metrics;


import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.benchmark.MultipleAlignment;
import org.biojava.bio.structure.align.model.AFPChain;

/**
 * Print the full length of one of the proteins (as opposed to the length
 * of the alignment, given by {@link AlignmentLengthMetric}
 * @author Spencer Bliven
 */
public class ProteinLengthMetric extends Metric {

	private int proteinIndex;

	public ProteinLengthMetric(int proteinIndex) {
		if(proteinIndex < 0 || proteinIndex > 1) {
			throw new IllegalArgumentException("proteinIndex must be 0 or 1");
		}

		this.proteinIndex = proteinIndex;
	}

	@Override
	public double calculate(MultipleAlignment reference, AFPChain align, Atom[] ca1, Atom[] ca2, Map<String, Object> metaData) {
		switch(proteinIndex) {
		case 0:
			return ca1.length;
		case 1:
			return ca2.length;
		default:
			//Never reached
			throw new IllegalStateException("Illegal Protein Index!");
		}
	}

	@Override
	public String getName() {
		return "Prot"+(proteinIndex)+"_len";
	}

	@Override
	public String format(double result) {
		return Integer.toString((int)result);
	}
}
