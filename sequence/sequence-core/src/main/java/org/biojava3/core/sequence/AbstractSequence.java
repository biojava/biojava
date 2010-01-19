package org.biojava3.core.sequence;

import java.util.Iterator;
import java.util.List;

import org.biojava3.core.sequence.impl.ArrayListSequenceBackingStore;

public abstract class AbstractSequence<C extends Compound> implements Sequence<C> {
	
	private final SequenceBackingStore<C> backingStore;
	
	public AbstractSequence(CharSequence seqString) {
		this(seqString, new ArrayListSequenceBackingStore<C>());
	}
	
	public AbstractSequence(CharSequence seqString, SequenceBackingStore<C> backingStore) {
		backingStore.setCompoundSet(this.getCompoundSet());
		backingStore.setContents(seqString);
		this.backingStore = backingStore;
	}
	
	public AbstractSequence(SequenceLazyLoader<C> lazyLoader) {
		this.backingStore = lazyLoader;
	}
	
	public abstract CompoundSet<C> getCompoundSet();
	
	public CharSequence getAsCharSequence() {
		return backingStore.getAsCharSequence();
	}

	public List<C> getAsList() {
		return backingStore.getAsList();
	}

	public C getCompoundAt(int position) {
		return backingStore.getCompoundAt(position);
	}

	public int getIndexOf(C compound) {
		return backingStore.getIndexOf(compound);
	}

	public int getLastIndexOf(C compound) {
		return backingStore.getLastIndexOf(compound);
	}

	public int getLength() {
		return backingStore.getLength();
	}

	public SequenceView<C> getSubSequence(final int start, final int end) {
		return new AbstractSequenceView<C>() {

			public int getEnd() {
				return end;
			}

			public int getStart() {
				return start;
			}

			public Sequence<C> getViewedSequence() {
				return AbstractSequence.this;
			}
		};
	}

	public Iterator<C> iterator() {
		return backingStore.iterator();
	}

}
