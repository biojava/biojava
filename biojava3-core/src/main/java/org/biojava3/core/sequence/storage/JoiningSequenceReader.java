package org.biojava3.core.sequence.storage;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.ProxySequenceReader;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.template.SequenceMixin;
import org.biojava3.core.sequence.template.SequenceView;

/**
 * This reader actually proxies onto multiple types of sequence in order
 * to allow a number of sequence objects to act as if they are one sequence.
 * The code takes in any number of sequences, records the minimum and maximum
 * bounds each sequence covers with respect to 1 position indexing and then
 * binary searches these when a position is requested. Because of this
 * 0 length Sequences are excluded during construction.
 *
 * Performance is not as good as if you are using a flat sequence however the
 * speed of lookup is more than adaquate for most situations. Using the iterator
 * gives the best performance as this does not rely on the binary search
 * mechanism instead iterating through each sequence in turn.
 *
 * @author ayates
 * @param <C> Tyoe of compound to hold
 */
public class JoiningSequenceReader<C extends Compound> implements ProxySequenceReader<C> {

    /**
     * Internal implementation flag and controls how we look for the right
     * sub-sequence
     */
    private static final boolean BINARY_SEARCH = true;
    private final List<Sequence<C>> sequences;
    private final CompoundSet<C> compoundSet;
    private int[] maxSequenceIndex;
    private int[] minSequenceIndex;

    /**
     * Allows creation of the store from Vargs Sequence<C> objects. CompoundSet
     * defaults to the first element of the array (assuming there are elements
     * available during construction otherwise we will throw an illegal
     * state exception).
     */
    public JoiningSequenceReader(Sequence<C>... sequences) {
        this(Arrays.asList(sequences));
    }

    /**
     * Allows creation of the store from List<Sequence<C>>. CompoundSet
     * defaults to the first element of the List (assuming there are elements
     * available during construction otherwise we will throw an illegal
     * state exception).
     */
    public JoiningSequenceReader(List<Sequence<C>> sequences) {
        this.sequences = grepSequences(sequences);
        this.compoundSet = grepCompoundSet();
    }

    public JoiningSequenceReader(CompoundSet<C> compoundSet, Sequence<C>... sequences) {
        this(compoundSet, Arrays.asList(sequences));
    }

    public JoiningSequenceReader(CompoundSet<C> compoundSet, List<Sequence<C>> sequences) {
        this.sequences = grepSequences(sequences);
        this.compoundSet = compoundSet;
    }

    private List<Sequence<C>> grepSequences(List<Sequence<C>> sequences) {
        List<Sequence<C>> seqs = new ArrayList<Sequence<C>>();
        for (Sequence<C> s : sequences) {
            if (s.getLength() != 0) {
                seqs.add(s);
            }
        }
        return seqs;
    }

    private CompoundSet<C> grepCompoundSet() {
        if (sequences.isEmpty()) {
            throw new IllegalStateException("Cannot get a CompoundSet because we have no sequences. Set during construction");
        }
        return sequences.get(0).getCompoundSet();
    }

    
    public C getCompoundAt(int position) {
        int sequenceIndex = getSequenceIndex(position);
        Sequence<C> sequence = sequences.get(sequenceIndex);
        int indexInSequence = (position - getMinSequenceIndex()[sequenceIndex]) + 1;
        return sequence.getCompoundAt(indexInSequence);
    }

    
    public CompoundSet<C> getCompoundSet() {
        return compoundSet;
    }

    
    public int getLength() {
        int[] maxSeqIndex = getMaxSequenceIndex();
        if (maxSeqIndex.length == 0) {
            return 0;
        }
        return maxSeqIndex[maxSeqIndex.length - 1];
    }

    /**
     * Returns which Sequence holds the position queried for
     */
    private int getSequenceIndex(int position) {
        if (BINARY_SEARCH) {
            return binarySearch(position);
        } else {
            return linearSearch(position);
        }
    }

    private int[] getMinSequenceIndex() {
        if (minSequenceIndex == null) {
            initSeqIndexes();
        }
        return minSequenceIndex;
    }

    private int[] getMaxSequenceIndex() {
        if (maxSequenceIndex == null) {
            initSeqIndexes();
        }
        return maxSequenceIndex;
    }

    private void initSeqIndexes() {
        minSequenceIndex = new int[sequences.size()];
        maxSequenceIndex = new int[sequences.size()];
        int currentMaxIndex = 0;
        int currentMinIndex = 1;
        int lastLength = 0;
        for (int i = 0; i < sequences.size(); i++) {
            currentMinIndex += lastLength;
            currentMaxIndex += sequences.get(i).getLength();
            minSequenceIndex[i] = currentMinIndex;
            maxSequenceIndex[i] = currentMaxIndex;
            lastLength = sequences.get(i).getLength();
        }
    }

    /**
     * Scans through the sequence index arrays in linear time. Not the best
     * performance but easier to code
     */
    private int linearSearch(int position) {
        int[] minSeqIndex = getMinSequenceIndex();
        int[] maxSeqIndex = getMaxSequenceIndex();
        int length = minSeqIndex.length;
        for (int i = 0; i < length; i++) {
            if (position >= minSeqIndex[i] && position <= maxSeqIndex[i]) {
                return i;
            }
        }
        throw new IndexOutOfBoundsException("Given position " + position + " does not map into this Sequence");
    }

    /**
     * Scans through the sequence index arrays using binary search
     */
    private int binarySearch(int position) {
        int[] minSeqIndex = getMinSequenceIndex();
        int[] maxSeqIndex = getMaxSequenceIndex();

        int low = 0;
        int high = minSeqIndex.length - 1;
        while (low <= high) {
            //Go to the mid point in the array
            int mid = (low + high) >>> 1;

            //Get the max position represented by this Sequence
            int midMinPosition = minSeqIndex[mid];
            int midMaxPosition = maxSeqIndex[mid];

            //if current position is greater than the current bounds then
            //increase search space
            if (midMinPosition < position && midMaxPosition < position) {
                low = mid + 1;
            } //if current position is less than current bounds then decrease
            //search space
            else if (midMinPosition > position && midMaxPosition > position) {
                high = mid - 1;
            } else {
                return mid;
            }
        }
        throw new IndexOutOfBoundsException("Given position " + position + " does not map into this Sequence");
    }

    /**
     * Iterator implementation which attempts to move through the 2D structure
     * attempting to skip onto the next sequence as & when it is asked to
     */
    
    public Iterator<C> iterator() {
        final List<Sequence<C>> localSequences = sequences;
        return new Iterator<C>() {

            private Iterator<C> currentSequenceIterator = null;
            private int currentPosition = 0;

            
            public boolean hasNext() {
                //If the current iterator is null then see if the Sequences object has anything
                if (currentSequenceIterator == null) {
                    return !localSequences.isEmpty();
                }

                //See if we had any compounds
                boolean hasNext = currentSequenceIterator.hasNext();
                if (!hasNext) {
                    hasNext = currentPosition < sequences.size();
                }
                return hasNext;
            }

            
            public C next() {
                if (currentSequenceIterator == null) {
                    if (localSequences.isEmpty()) {
                        throw new NoSuchElementException("No sequences to iterate over; make sure you call hasNext() before next()");
                    }
                    currentSequenceIterator = localSequences.get(currentPosition).iterator();
                    currentPosition++;
                }
                if (!currentSequenceIterator.hasNext()) {
                    currentSequenceIterator = localSequences.get(currentPosition).iterator();
                    currentPosition++;
                }
                return currentSequenceIterator.next();
            }

            
            public void remove() throws UnsupportedOperationException {
                throw new UnsupportedOperationException("Cannot remove from this Sequence");
            }
        };
    }

    
    public void setCompoundSet(CompoundSet<C> compoundSet) {
        throw new UnsupportedOperationException();
    }

    
    public void setContents(String sequence) {
        throw new UnsupportedOperationException();
    }

    
    public int countCompounds(C... compounds) {
        return SequenceMixin.countCompounds(this, compounds);
    }

    
    public AccessionID getAccession() throws UnsupportedOperationException {
        throw new UnsupportedOperationException();
    }

    
    public List<C> getAsList() {
        return SequenceMixin.toList(this);
    }

    
    public int getIndexOf(C compound) {
        return SequenceMixin.indexOf(this, compound);
    }

    
    public int getLastIndexOf(C compound) {
        return SequenceMixin.lastIndexOf(this, compound);
    }

    
    public String getSequenceAsString() {
        return SequenceMixin.toStringBuilder(this).toString();
    }

    
    public SequenceView<C> getSubSequence(Integer start, Integer end) {
        return SequenceMixin.createSubSequence(this, start, end);
    }

    @Override
    public SequenceView<C> getInverse() {
        return SequenceMixin.inverse(this);
    }
}
