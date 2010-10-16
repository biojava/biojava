package org.biojava3.core.sequence.views;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

/**
 * A sliding window view of a sequence which does not implement any
 * interfaces like {@link Sequence} because they do not fit how this works.
 * For each index requested we return a List of compounds back.
 *
 * If you perform a view on a Sequence whose length is not a multiple of the
 * window the final window will be omitted i.e. if we have the sequence AGCGG
 * and a window of 3 then you will only see AGC since GG exceeds the calculated
 * length of this sequence.
 *
 * Because this does not implement a Sequence interface we do not recommend
 * passing this class around. If you need to represent a windowed sequence
 * as a real Sequence then translate it into a new Compound
 *
 * @author ayates
 *
 * @param <C> The type of compound we return from a window
 */
public class WindowedSequence<C extends Compound> implements Iterable<List<C>> {

    private final Sequence<C> sequence;
    private final int windowSize;

    public WindowedSequence(Sequence<C> sequence, int windowSize) {
        this.sequence = sequence;
        this.windowSize = windowSize;
    }

    /**
     * Access the current window size
     */
    public int getWindowSize() {
        return windowSize;
    }

    /**
     * Access the sequence which backs this window
     */
    public Sequence<C> getBackingSequence() {
        return sequence;
    }

    /**
     * Calculates start index according to the equation start = ( (index-1) -
     * windowSize) +1
     */
    protected int toStartIndex(int index) {
        return ((index - 1) * getWindowSize()) + 1;
    }

    /**
     * Returns the size of the windowed sequence which is the length by the
     * window size. Trailing Compounds are omitted.
     */
    public int getLength() {
        return getBackingSequence().getLength() / getWindowSize();
    }

    /**
     * For a given position into the windowed view this will return those
     * compounds we can see in the window. i.e. in the sequence AGGCCT requesting
     * index 1 returns AGG and requesting index 2 return CCT.
     *
     * @param index Windowed index position
     * @return The List of compounds
     */
    public List<C> getCompounds(int index) {
        int start = toStartIndex(index);
        int window = getWindowSize();
        Sequence<C> backing = getBackingSequence();
        List<C> output = new ArrayList<C>(window);
        for (int i = 0; i < window; i++) {
            output.add(backing.getCompoundAt(start + i));
        }
        return output;
    }

    /**
     * Iterator of all List of compounds available in a windowed sequence.
     */
    @Override
    public Iterator<List<C>> iterator() {
        return new WindowedIterator<C>(this);
    }

    /**
     * Used rather than an inline iterator which does its own iteration
     * not delegating back to the original class due to speed
     */
    private static class WindowedIterator<C extends Compound> implements Iterator<List<C>> {

        private final int end;
        private final int window;
        private int currentIndex = 1;
        private final Iterator<C> iter;

        public WindowedIterator(WindowedSequence<C> sequence) {
            this.window = sequence.getWindowSize();
            this.iter = sequence.getBackingSequence().iterator();
            this.end = sequence.getLength();
        }

        @Override
        public boolean hasNext() {
            return currentIndex <= end;
        }

        @Override
        public List<C> next() {
            List<C> results = new ArrayList<C>(window);
            for(int i = 0; i < window; i++) {
                results.add(iter.next());
            }
            currentIndex++;
            return results;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("Cannot remove from a Windowed view");
        }
    }
}
