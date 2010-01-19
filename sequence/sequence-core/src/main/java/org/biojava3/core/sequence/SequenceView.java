package org.biojava3.core.sequence;

public interface SequenceView<C extends Compound> extends Sequence<C> {
	
	public Sequence<C> getViewedSequence();
	
	/**
	 * 1-indexed, inclusive.
	 */
	public int getStart();
	
	/**
	 * 1-indexed, inclusive.
	 */
	public int getEnd();
}
