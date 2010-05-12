package org.biojava3.core.sequence.template;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.storage.ArrayListSequenceReader;
import org.biojava3.core.sequence.views.ReversedSequenceView;
import org.biojava3.core.util.CRC64Checksum;

import com.google.common.base.Function;
import com.google.common.collect.MapMaker;

/**
 * Provides a set of static methods to be used as static imports when needed
 * across multiple Sequence implementations but inheritance gets in the way.
 *
 * It also provides a place to put utility methods whose application can
 * be to a single class of Sequence e.g. {@link NucleotideCompound}
 * {@link Sequence}; or to any Sequence e.g. looking for the
 * {@link #getComposition(Sequence)} or {@link #getDistribution(Sequence)}
 * for any type of Sequence.
 *
 * All of these methods assume that you can use the {@link Iterable} interface
 * offered by the implementations of {@link Sequence} to provide all the
 * compounds that implementation allows you to see. Since sequence should know
 * nothing about its backing stores (apart from calling out to it) this should
 * be true.
 *
 * @author ayates
 */
public class SequenceMixin {

    /**
     * For the given vargs of compounds this method counts the number of
     * times those compounds appear in the given sequence
     *
     * @param sequence The {@link Sequence} to perform the count on
     * @param compounds The compounds to look for
     * @param <C> The type of compound we are looking for
     * @return The number of times the given compounds appear in this Sequence
     */
    public static <C extends Compound> int countCompounds(
            Sequence<C> sequence, C... compounds) {
        int count = 0;
        Map<C, Integer> compositon = getComposition(sequence);
        for (C compound : compounds) {
            count = compositon.get(compound) + count;
        }
        return count;
    }

    /**
     * Returns the count of GC in the given sequence
     *
     * @param sequence The {@link NucleotideCompound} {@link Sequence} to perform
     * the GC analysis on
     * @return The number of GC compounds in the sequence
     */
    public static int countGC(Sequence<NucleotideCompound> sequence) {
        CompoundSet<NucleotideCompound> cs = sequence.getCompoundSet();
        NucleotideCompound G = cs.getCompoundForString("G");
        NucleotideCompound C = cs.getCompoundForString("C");
        NucleotideCompound g = cs.getCompoundForString("g");
        NucleotideCompound c = cs.getCompoundForString("c");
        return countCompounds(sequence, G, C, g, c);
    }

    /**
     * Returns the count of AT in the given sequence
     *
     * @param sequence The {@link NucleotideCompound} {@link Sequence} to perform
     * the AT analysis on
     * @return The number of AT compounds in the sequence
     */
    public static int countAT(Sequence<NucleotideCompound> sequence) {
        CompoundSet<NucleotideCompound> cs = sequence.getCompoundSet();
        NucleotideCompound A = cs.getCompoundForString("A");
        NucleotideCompound T = cs.getCompoundForString("T");
        NucleotideCompound a = cs.getCompoundForString("a");
        NucleotideCompound t = cs.getCompoundForString("t");
        return countCompounds(sequence, A, T, a, t);
    }

    /**
     * Analogous to {@link #getComposition(Sequence)} but returns the
     * distribution of that {@link Compound} over the given sequence.
     *
     * @param <C> The type of compound to look for
     * @param sequence The type of sequence to look over
     * @return Returns the decimal fraction of the compounds in the given
     * sequence. Any compound not in the Map will return a fraction of 0.
     */
    public static <C extends Compound> Map<C, Double> getDistribution(Sequence<C> sequence) {
        Map<C, Double> results = new MapMaker().makeComputingMap(new Function<C, Double>() {

            public Double apply(C compound) {
                return 0.0;
            }
        });
        Map<C, Integer> composition = getComposition(sequence);
        double length = (double) sequence.getLength();
        for (Map.Entry<C, Integer> entry : composition.entrySet()) {
            double dist = entry.getValue().doubleValue() / length;
            results.put(entry.getKey(), dist);
        }
        return results;
    }

    /**
     * Does a linear scan over the given Sequence and records the number of
     * times each base appears. The returned map will return 0 if a compound
     * is asked for and the Map has no record of it.
     *
     * @param <C> The type of compound to look for
     * @param sequence The type of sequence to look over
     * @return Counts for the instances of all compounds in the sequence
     */
    public static <C extends Compound> Map<C, Integer> getComposition(Sequence<C> sequence) {
        Map<C, Integer> results = new MapMaker().makeComputingMap(new Function<C, Integer>() {

            public Integer apply(C compound) {
                return 0;
            }
        });

        for (C currentCompound : sequence) {
            Integer currentInteger = results.get(currentCompound);
            currentInteger++;
            results.put(currentCompound, currentInteger);
        }
        return results;
    }

    /**
     * For the given Sequence this will return a {@link StringBuilder} object
     * filled with the results of {@link Compound#toString()}.
     */
    public static <C extends Compound> StringBuilder toStringBuilder(Sequence<C> sequence) {
        StringBuilder sb = new StringBuilder(sequence.getLength());
        for (C compound : sequence) {
            sb.append(compound.toString());
        }
        return sb;
    }

    /**
     * For the given {@link Sequence} this will return a {@link List} filled with
     * the Compounds of that {@link Sequence}.
     */
    public static <C extends Compound> List<C> toList(Sequence<C> sequence) {
        List<C> list = new ArrayList<C>(sequence.getLength());
        for (C compound : sequence) {
            list.add(compound);
        }
        return list;
    }

    /**
     * Performs a linear search of the given Sequence for the given compound.
     * Once we find the compound we return the position.
     */
    public static <C extends Compound> int indexOf(Sequence<C> sequence,
            C compound) {
        int index = 1;
        for (C currentCompound : sequence) {
            if (currentCompound.equals(compound)) {
                return index;
            }
            index++;
        }
        return 0;
    }

    /**
     * Performs a reversed linear search of the given Sequence by wrapping
     * it in a {@link ReversedSequenceView} and passing it into
     * {@link #indexOf(Sequence, Compound)}.
     */
    public static <C extends Compound> int lastIndexOf(Sequence<C> sequence,
            C compound) {
        return indexOf(new ReversedSequenceView<C>(sequence), compound);
    }

    /**
     * Creates a simple sequence iterator which moves through a sequence going
     * from 1 to the length of the Sequence. Modification of the Sequence is not
     * allowed.
     */
    public static <C extends Compound> Iterator<C> createIterator(
            Sequence<C> sequence) {
        return new SequenceIterator<C>(sequence);
    }

    /**
     * Creates a simple sub sequence view delimited by the given start and end.
     */
    public static <C extends Compound> SequenceView<C> createSubSequence(
            Sequence<C> sequence, int start, int end) {
        return new SubSequence<C>(sequence, start, end);
    }

    /**
     * Implements sequence shuffling by first materializing the given
     * {@link Sequence} into a {@link List}, applying
     * {@link Collections#shuffle(List)} and then returning the shuffled
     * elements in a new instance of {@link SequenceBackingStore} which behaves
     * as a {@link Sequence}.
     */
    public static <C extends Compound> Sequence<C> shuffle(Sequence<C> sequence) {
        List<C> compounds = sequence.getAsList();
        Collections.shuffle(compounds);
        return new ArrayListSequenceReader<C>(compounds,
                sequence.getCompoundSet());
      
    }

    /**
     * Performs a simple CRC64 checksum on any given sequence.
     */
    public static <C extends Compound> String checksum(Sequence<C> sequence) {
        CRC64Checksum checksum = new CRC64Checksum();
        for (C compound : sequence) {
            checksum.update(compound.getShortName());
        }
        return checksum.toString();
    }

    /**
     * A basic sequence iterator which iterates over the given Sequence by
     * biological index. This assumes your sequence supports random access
     * and performs well when doing these operations.
     *
     * @author ayates
     *
     * @param <C> Type of compound to return
     */
    private static class SequenceIterator<C extends Compound>
            implements Iterator<C> {

        private final Sequence<C> sequence;
        private final int length;
        private int currentPosition = 0;

        public SequenceIterator(Sequence<C> sequence) {
            this.sequence = sequence;
            this.length = sequence.getLength();
        }

        public boolean hasNext() {
            if (currentPosition < length) {
                currentPosition++;
                return true;
            }
            return false;
        }

        public C next() {
            return sequence.getCompoundAt(currentPosition);
        }

        public void remove() {
            throw new UnsupportedOperationException("Cannot remove() on a SequenceIterator");
        }
    }

    /**
     * Simple sub-sequence object which holds a start and a end value which
     * allows us to look at subsequences for the given Sequence. The sub-sequence
     * is then referred to with respect to that sequence e.g. creating a
     * sub sequence which starts at 5 and ends at 10; when querying for position
     * 1 in the sub sequence we get back 5, when querying for 2 we get back 6.
     *
     * @author ayates
     *
     * @param <C> Type of compound to hold
     */
    private static class SubSequence<C extends Compound>
            extends SequenceProxyView<C> {

        private final Sequence<C> viewedSequence;
        private final int start;
        private final int end;



 //           public SubSequence( Sequence<C> sequence){
 //       this(sequence,1,sequence.getLength());
 //   }

        public SubSequence(Sequence<C> viewedSequence, Integer start, Integer end) {
            this.viewedSequence = viewedSequence;
            this.start = start;
            this.end = end;
        }

        public int getEnd() {
            return end;
        }

        public int getStart() {
            return start;
        }

        public Sequence<C> getViewedSequence() {
            return viewedSequence;
        }

        @Override
        public AccessionID getAccession() {
            return getViewedSequence().getAccession();
        }
    }
}
