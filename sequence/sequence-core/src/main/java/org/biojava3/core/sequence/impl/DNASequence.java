package org.biojava3.core.sequence.impl;

import org.biojava3.core.sequence.AbstractSequence;
import org.biojava3.core.sequence.CompoundSet;
import org.biojava3.core.sequence.SequenceBackingStore;
import org.biojava3.core.sequence.SequenceLazyLoader;

public class DNASequence extends AbstractSequence<DNACompound> {

	private final CompoundSet<DNACompound> compoundSet = new DNACompoundSet();
	
	public DNASequence(CharSequence seqString,
			SequenceBackingStore<DNACompound> backingStore) {
		super(seqString, backingStore);
	}

	public DNASequence(CharSequence seqString) {
		super(seqString);
	}

	public DNASequence(SequenceLazyLoader<DNACompound> lazyLoader) {
		super(lazyLoader);
	}

	public CompoundSet<DNACompound> getCompoundSet() {
		return compoundSet;
	}

}
