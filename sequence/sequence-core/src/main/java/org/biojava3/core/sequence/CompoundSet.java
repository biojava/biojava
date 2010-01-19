package org.biojava3.core.sequence;

public interface CompoundSet<C extends Compound> {
	public int getMaxSingleCompoundCharSequenceLength();
	
	/**
	 * Return null if not recognised. Throw IllegalArgumentException if string
	 * is longer than maximum allowed by {@link #getCharSequenceForCompound(Compound)}.
	 */
	public C getCompoundForCharSequence(CharSequence string);
	
	public CharSequence getCharSequenceForCompound(C compound);
}
