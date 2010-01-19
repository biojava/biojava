package org.biojava3.core.sequence;

import java.util.List;

public interface Sequence<C extends Compound> extends Iterable<C> {
	public int getLength();
	
	public C getCompoundAt(int position);
	
	public int getIndexOf(C compound);
	
	public int getLastIndexOf(C compound);
	
	public CharSequence getAsCharSequence();
	
	public List<C> getAsList();

	public SequenceView<C> getSubSequence(int start, int end);
}
