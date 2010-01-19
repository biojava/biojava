package org.biojava3.core.sequence;

public interface SequenceBackingStore<C extends Compound> extends Sequence<C> {
	
	public void setCompoundSet(CompoundSet<C> compoundSet);
	
	public void setContents(CharSequence charSequence);
}
