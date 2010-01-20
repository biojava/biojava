package org.biojava3.core.sequence.template;

public interface CompoundSet<C extends Compound> {
	public int getMaxSingleCompoundStringLength();
	
	/**
	 * Return null if not recognised. Throw IllegalArgumentException if string
	 * is longer than maximum allowed by {@link #getStringForCompound(Compound)}.
	 */
	public C getCompoundForString(String string);
	
	public String getStringForCompound(C compound);
}
