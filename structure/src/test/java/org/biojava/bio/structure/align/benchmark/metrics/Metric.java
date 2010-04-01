/**
 * 
 */
package org.biojava.bio.structure.align.benchmark.metrics;


import org.biojava.bio.structure.align.benchmark.MultipleAlignment;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.Atom;

/**
 * Represents some measure of comparison between a calculated alignment
 * (eg the result of some {@link StructureAlignment#align(Atom[],Atom[])} call)
 * with a reference {@link Alignment}.
 * <p>
 * Metrics should be written as singletons. They may make use of instance data
 * (eg parameters, etc), but the calculate() and toString() methods should not
 * modify the instance data, such that multiple calls to {@link #calculate(MultipleAlignment, AFPChain, Atom[], Atom[])}
 * with the same parameters should return the same result.
 * 
 * @author Spencer Bliven
 */
public abstract class Metric {

	/**
	 * Calculate the Metric for a particular alignment & reference.
	 * 
	 * @param reference The reference alignment against which to judge align
	 * @param align The test alignment, to compare against the reference
	 * @param ca1 TODO
	 * @param ca2 TODO
	 * @return
	 */
	public abstract double calculate(MultipleAlignment reference, AFPChain align, Atom[] ca1, Atom[] ca2);
	
	/**
	 * Calculates the Metric for the given parameters and returns the string
	 * representation of that object. Equivalent to calling
	 * <code>this.toString(calculate(reference, align, ca1, ca2))</code>
	 * @param reference The reference alignment against which to judge align.
	 * @param align The test alignment, to compare against the reference
	 * @param ca1 The CA atoms of the first 
	 * @param ca2 TODO
	 * @return A string representation of the Object returned by calculate()
	 */
	public String format(MultipleAlignment reference, AFPChain align, Atom[] ca1, Atom[] ca2) {
		return this.format(calculate(reference, align, ca1, ca2));
	}
	
	/**
	 * To facilitate outputting, this method converts the Object returned by 
	 * {@link #calculate(Alignment, AFPChain, Atom[], Atom[]) this.calculate(goldStandard, align)}
	 * into a String. The default implementation is just to call the
	 * {@link Object#toString()} method, but in some metrics a shorter
	 * representation may be more appropriate for tab-delimited output.
	 * <p>
	 * The output of this method should not contain tab characters.
	 * @param result An object returned by a previous call to calculate()
	 * @return The string representation of that object
	 */
	public String format(double result) {
		return Double.toString(result);
	}
	
	/**
	 * @return A one-word description of this Metric, suitable for use as a header in an
	 * output file. Should not contain whitespace.
	 */
	public abstract String getName();
	
	/**
	 * Default toString; same as getName().
	 */
	public String toString() { return getName(); }
}
