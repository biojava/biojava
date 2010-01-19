package org.biojava3.core.sequence.impl;

import org.biojava3.core.sequence.AbstractSequence;
import org.biojava3.core.sequence.CompoundSet;
import org.biojava3.core.sequence.SequenceBackingStore;
import org.biojava3.core.sequence.SequenceLazyLoader;

public class GenericSequence extends AbstractSequence<GenericCompound> {

	private final CompoundSet<GenericCompound> compoundSet = new GenericCompoundSet();
	
	public GenericSequence(CharSequence seqString,
			SequenceBackingStore<GenericCompound> backingStore) {
		super(seqString, backingStore);
	}

	public GenericSequence(CharSequence seqString) {
		super(seqString);
	}

	public GenericSequence(SequenceLazyLoader<GenericCompound> lazyLoader) {
		super(lazyLoader);
	}

	public CompoundSet<GenericCompound> getCompoundSet() {
		return compoundSet;
	}

}
