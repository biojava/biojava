package org.biojava3.core.sequence;

public class CompoundNotFoundError extends Error {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public CompoundNotFoundError(CharSequence compoundStr) {
		super("Compound not found for: "+compoundStr);
	}
}
