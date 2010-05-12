package org.biojava3.core.sequence.storage;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.template.SequenceMixin;
import org.biojava3.core.sequence.template.ProxySequenceReader;
import org.biojava3.core.sequence.template.SequenceView;

/**
 * An implementation of the popular 2bit encoding. 2bit is the process
 * of representing ATGC as a series of states in 2 bits. For example
 * we can encode the base T as 00 (2 bits in the off state). This means
 * that for every byte we can store 4 bases (since there are 8 bits in
 * a byte). The actual storage is dependent on subclasses of
 * {@link TwoBitArrayWorker} but normally we work using a backing data
 * structure of int[] so we actually store 16 bases for every
 * int value.
 *
 * This should be more memory efficient than using the
 * {@link ArrayListSequenceBackingStore} but there are a number of issues
 * with using this format.
 *
 * <ul>
 * <li>We can only support 4 {@link Compound}s so no N</li>
 * <li>For real savings you must read the sequence in using your own
 * Reader and a {@link TwoBitArrayWorker} instance</li>
 * </ul>
 *
 * @author ayates
 *
 * @param <C> Type of compound; must extend {@link NucleotideCompound}
 */
public class TwoBitSequenceBackingStore<C extends NucleotideCompound> implements ProxySequenceReader<C> {

    private final AccessionID accession;
    private final TwoBitArrayWorker<C> worker;

    public TwoBitSequenceBackingStore(Sequence<C> sequence) {
        this.accession = sequence.getAccession();
        worker = new TwoBitArrayWorker<C>(sequence.getCompoundSet(), sequence.getLength());
        worker.populate(sequence);
    }

    public TwoBitSequenceBackingStore(String sequence, CompoundSet<C> compoundSet) {
        this(sequence, compoundSet, new AccessionID("Unknown"));
    }

    public TwoBitSequenceBackingStore(String sequence, CompoundSet<C> compoundSet, AccessionID accession) {
        worker = new TwoBitArrayWorker<C>(compoundSet, sequence.length());
        this.accession = accession;
        worker.populate(sequence);
    }

    public TwoBitSequenceBackingStore(TwoBitArrayWorker<C> worker, AccessionID accession) {
        this.accession = accession;
        this.worker = worker;
    }

    /**
     * Class is immutable & so this is unsupported
     */
    public void setCompoundSet(CompoundSet<C> compoundSet) {
        throw new UnsupportedOperationException("Cannot reset the CompoundSet; object is immutable");
    }

    /**
     * Class is immutable & so this is unsupported
     */
    public void setContents(String sequence) {
        throw new UnsupportedOperationException(this.getClass().getSimpleName() + " is an immutable data structure; cannot reset contents");
    }

    /**
     * Counts the number of times a compound appears in this sequence store
     */
    public int countCompounds(C... compounds) {
        return SequenceMixin.countCompounds(this, compounds);
    }

    public AccessionID getAccession() {
        return accession;
    }

    /**
     * Returns this Sequence store as a List
     */
    public List<C> getAsList() {
        return SequenceMixin.toList(this);
    }

    /**
     * Returns the compound at the specified biological index
     */
    public C getCompoundAt(int position) {
        return worker.getCompoundAt(position);
    }

    /**
     * Returns the compound set backing this store
     */
    public CompoundSet<C> getCompoundSet() {
        return worker.getCompoundSet();
    }

    /**
     * Returns the first occurrence of the given compound in this store; performs
     * a linear search
     */
    public int getIndexOf(C compound) {
        return SequenceMixin.indexOf(this, compound);
    }

    /**
     * Returns the last occurrence of the given compound in this store; performs
     * a linear search
     */
    public int getLastIndexOf(C compound) {
        return SequenceMixin.lastIndexOf(this, compound);
    }

    /**
     * Returns the length of the sequence
     */
    public int getLength() {
        return worker.getLength();
    }

    /**
     * Returns the sequence as a String
     */
    public String getSequenceAsString() {
        return SequenceMixin.toStringBuilder(this).toString();
    }

    /**
     * Returns a sub sequence view
     */
    public SequenceView<C> getSubSequence(final int start, final int end) {
        return SequenceMixin.createSubSequence(this, start, end);
    }

    /**
     * Provides basic iterable access to this class
     */
    public Iterator<C> iterator() {
        return SequenceMixin.createIterator(this);
    }

    @Override
    public String getSequenceAsString(Integer start, Integer end, Strand strand) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public SequenceView<C> getSubSequence(Integer start, Integer end) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * The logic of working with 2Bit has been separated out into this class
     * to help developers create 2Bit data structures without having to
     * put the code into an intermediate format and to also use the format
     * without the need to copy this code.
     *
     * This version works with a
     *
     * This class behaves just like a {@link Sequence} without the interface
     *
     * @author ayates
     *
     * @param <C> The {@link Compound} to use; should subclass
     * {@link NucleotideCompound}
     */
    public static class TwoBitArrayWorker<C extends NucleotideCompound> {

        private final CompoundSet<C> compoundSet;
        private final int length;
        private final int[] sequence;
        private transient List<C> indexToCompoundsLookup = null;
        private transient Map<C, Integer> compoundsToIndexLookup = null;
        public static final int BYTES_PER_DATATYPE = 32;
        /**
         * Masking value used for extracting the right most 2 bits from a byte
         */
        private final static byte MASK = (byte) ((int) Math.pow(2, 0) | (int) Math.pow(2, 1));

        public TwoBitArrayWorker(CompoundSet<C> compoundSet, int length) {
            this.compoundSet = compoundSet;
            this.length = length;
            this.sequence = new int[seqArraySize(length)];
        }

        public TwoBitArrayWorker(CompoundSet<C> compoundSet, int[] sequence) {
            this.compoundSet = compoundSet;
            this.sequence = sequence;
            this.length = sequence.length;
        }

        public static int seqArraySize(int length) {
            return (int) Math.ceil((double) length / (double) BYTES_PER_DATATYPE);
        }

        /**
         * Loops through the Compounds in a Sequence and passes them onto
         * {@link #setCompoundAt(Compound, int)}
         */
        public void populate(Sequence<C> sequence) {
            int position = 1;
            for (C c : sequence) {
                setCompoundAt(c, position++);
            }
        }

        /**
         * Loops through the chars in a String and passes them onto
         * {@link #setCompoundAt(char, int)}
         */
        public void populate(String sequence) {
            for (int index = 0; index < getLength(); index++) {
                setCompoundAt(sequence.charAt(index), index + 1);
            }
        }

        /**
         * Converts from char to Compound and sets it at the given biological index
         */
        public void setCompoundAt(char base, int position) {
            C compound = getCompoundSet().getCompoundForString(Character.toString(base));
            setCompoundAt(compound, position);
        }

        /**
         * Sets the compound at the specified biological index
         */
        public void setCompoundAt(C compound, int position) {
            int arrayIndex = biologicalIndexToArrayIndex(position);
            int currentInt = sequence[arrayIndex];
            int shiftBy = shiftBy(position);
            Integer integerValue = getCompoundsToIndexLookup().get(compound);

            //If we got nothing then throw an error as it's wrong
            if (integerValue == null) {
                processUnknownCompound(compound, position);
            }

            int shiftedValue = integerValue << shiftBy;

            sequence[arrayIndex] = currentInt | shiftedValue;
        }

        /**
         * Returns the compound at the specified biological index
         */
        public C getCompoundAt(int position) {
            //Avoids asking for something which is not encoded by a bit-pair
            if (position > getLength()) {
                throw new IllegalArgumentException(position + " is greater than length. Cannot access this position");
            }
            //Just stops us from using 0 indexing
            if (position < 1) {
                throw new IllegalArgumentException(position + " is less than 1; you must use biological indexing (indexing from 1)");
            }

            int arrayIndex = biologicalIndexToArrayIndex(position);
            int currentByte = sequence[arrayIndex];
            int shiftBy = shiftBy(position);
            int shifted = (int) (currentByte >>> shiftBy);
            int masked = (int) (shifted & MASK);

            if (masked > 3) {
                throw new IllegalStateException("Got a masked value of " + masked + "; do not understand values greater than 3");
            }
            return getIndexToCompoundsLookup().get(masked);
        }

        /**
         * Since 2bit only supports 4 bases it is more than likely when processing
         * sequence you will encounter a base which is not one of the four (such
         * as an N). Currently this method will throw an exception detailing the
         * base it cannot process.
         *
         * You can override this to convert the unknown base into one you can
         * process or store locations of unknown bases for a level of post processing
         * in your subclass.
         *
         * @param compound Compound process
         * @return Byte representation of the compound
         * @throws IllegalStateException Done whenever this method is invoked
         */
        protected byte processUnknownCompound(C compound, int position) throws IllegalStateException {
            throw new IllegalStateException("Do not know how to translate the compound " + compound + " to a 2bit representation");
        }

        /**
         * Returns a list of compounds the index position of which is used
         * to translate from the byte representation into a compound.
         */
        @SuppressWarnings("serial")
        protected List<C> getIndexToCompoundsLookup() {
            if (indexToCompoundsLookup == null) {
                final CompoundSet<C> cs = getCompoundSet();
                indexToCompoundsLookup = new ArrayList<C>() {

                    {
                        add(cs.getCompoundForString("T"));
                        add(cs.getCompoundForString("C"));
                        add(cs.getCompoundForString("A"));
                        add(cs.getCompoundForString("G"));
                    }
                };
            }
            return indexToCompoundsLookup;
        }

        /**
         * Returns a map which converts from compound to byte
         */
        @SuppressWarnings("serial")
        protected Map<C, Integer> getCompoundsToIndexLookup() {
            if (compoundsToIndexLookup == null) {
                final CompoundSet<C> cs = getCompoundSet();
                compoundsToIndexLookup = new HashMap<C, Integer>() {

                    {
                        put(cs.getCompoundForString("T"), 0);
                        put(cs.getCompoundForString("C"), 1);
                        put(cs.getCompoundForString("A"), 2);
                        put(cs.getCompoundForString("G"), 3);
                        put(cs.getCompoundForString("t"), 0);
                        put(cs.getCompoundForString("c"), 1);
                        put(cs.getCompoundForString("a"), 2);
                        put(cs.getCompoundForString("g"), 3);
                    }
                };
            }
            return compoundsToIndexLookup;
        }

        /**
         * Array position is biological - 1 / 4
         * e.g. for position 4 we store this in array index 0 & BYTES_PER_DATATYPE = 4
         * (4-1) = 3
         * 3/4 = 0.75 (floor to 0)
         */
        private int biologicalIndexToArrayIndex(int index) {
            return ((index - 1) / BYTES_PER_DATATYPE);
        }

        /**
         * Convert from bio to 0 index, remainder & then multiply by 2
         * e.g. position 7 is the 3rd position in a byte but a shift of 4 therefore
         * 7 - 1 = 6
         * 6 % 4 = 2
         * 2 * 2 = 4
         */
        private byte shiftBy(int index) {
            return (byte) (((index - 1) % BYTES_PER_DATATYPE) * 2);
        }

        /**
         * Returns the compound set backing this store
         */
        public CompoundSet<C> getCompoundSet() {
            return compoundSet;
        }

        public int getLength() {
            return length;
        }
    }
}
