package org.biojava.bio.structure.align.benchmark.metrics;


import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.benchmark.MultipleAlignment;
import org.biojava.bio.structure.align.model.AFPChain;

public abstract class LengthMetric {
	/**
	 * Calculates the length of the reference alignment
	 * @author Spencer Bliven
	 *
	 */
	public static class Reference extends Metric {

		@Override
		public double calculate(MultipleAlignment reference, AFPChain align, Atom[] ca1, Atom[] ca2) {
			return reference.size();
		}

		@Override
		public String getName() {
			return "Ref_len";
		}

		@Override
		public String format(double result) {
			return Integer.toString((int)result);
		}
	}
	
	/**
	 * Calculates the length of the test alignment
	 * @author Spencer Bliven
	 *
	 */
	public static class Alignment extends Metric {

		@Override
		public double calculate(MultipleAlignment reference, AFPChain align, Atom[] ca1, Atom[] ca2) {
			return (double)align.getAlnLength();
		}

		@Override
		public String getName() {
			return "Aln_len";
		}
		
		@Override
		public String format(double result) {
			return Integer.toString((int)result);
		}
	}
}
