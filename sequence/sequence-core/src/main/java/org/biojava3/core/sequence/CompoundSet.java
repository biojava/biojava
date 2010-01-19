package org.biojava3.core.sequence;

public interface CompoundSet<C extends Compound> {
	public int getMaxSingleCompoundCharSequenceLength();
	
	public C getCompoundForCharSequence(CharSequence string);
	
	public CharSequence getCharSequenceForCompound(C compound);
}
