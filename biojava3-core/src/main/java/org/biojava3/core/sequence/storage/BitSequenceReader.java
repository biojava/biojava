package org.biojava3.core.sequence.storage;

import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.SequenceMixin;
import org.biojava3.core.sequence.template.ProxySequenceReader;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.template.SequenceView;
import org.biojava3.core.util.Equals;
import org.biojava3.core.util.Hashcoder;

/**
 * An implementation of the popular bit encodings. This class provides the
 * Sequence view over what is actually carried out in the {@link BitArrayWorker}
 * instances. These are the objects that carry out array storage as well as
 * indexing into those arrays. New bit encodings can be written by extending
 * this class and a worker class. There are a number of issues with this
 * type of storage engine:
 *
 * <ul>
 * <li>We can only support a finite number of {@link Compound}s; 2 bit allows no N compounds</li>
 * <li>For real savings you must read the sequence in using your own
 * Reader and a {@link BitArrayWorker} instance</li>
 * </ul>
 *
 * @author ayates
 *
 * @param <C> Type of compound; must extend {@link NucleotideCompound}
 */
public class BitSequenceReader<C extends Compound> implements ProxySequenceReader<C> {

    private final AccessionID accession;
    private final BitArrayWorker<C> worker;

    /**
     * Instance which allows you to supply a different @{BitArrayWorker}
     * object.
     */
    public BitSequenceReader(BitArrayWorker<C> worker, AccessionID accession) {
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
        throw new UnsupportedOperationException(getClass().getSimpleName() + " is an immutable data structure; cannot reset contents");
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
    @Override
    public Iterator<C> iterator() {
        return SequenceMixin.createIterator(this);
    }
    
    public SequenceView<C> getSubSequence(Integer start, Integer end) {
        return getSubSequence((int) start, (int) end);
    }

    @Override
    public SequenceView<C> getInverse() {
        return SequenceMixin.inverse(this);
    }

    @Override
    public int hashCode() {
        int s = Hashcoder.SEED;
        s = Hashcoder.hash(s, accession);
        s = Hashcoder.hash(s, worker);
        return s;
    }

    @Override
    public boolean equals(Object o) {
        if(Equals.classEqual(this, o)) {
            @SuppressWarnings("unchecked")
            BitSequenceReader<C> that = (BitSequenceReader<C>)o;
            return  Equals.equal(this.accession, that.accession) &&
                    Equals.equal(this.worker, that.worker);
        }
        return false;
    }

    /**
     * The logic of working with a bit has been separated out into this class
     * to help developers create the bit data structures without having to
     * put the code into an intermediate format and to also use the format
     * without the need to copy this code.
     *
     * This class behaves just like a {@link Sequence} without the interface
     *
     * @author ayates
     *
     * @param <C> The {@link Compound} to use
     */
    public static abstract class BitArrayWorker<C extends Compound> {

        private final CompoundSet<C> compoundSet;
        private final int length;
        private final int[] sequence;
        private transient List<C> indexToCompoundsLookup = null;
        private transient Map<C, Integer> compoundsToIndexLookup = null;
        public static final int BYTES_PER_INT = 32;

        private volatile Integer hashcode = null;

        public BitArrayWorker(Sequence<C> sequence) {
            this(sequence.getCompoundSet(), sequence.getLength());
            populate(sequence);
        }

        public BitArrayWorker(String sequence, CompoundSet<C> compoundSet) {
            this(compoundSet, sequence.length());
            populate(sequence);
        }

        public BitArrayWorker(CompoundSet<C> compoundSet, int length) {
            this.compoundSet = compoundSet;
            this.length = length;
            this.sequence = new int[seqArraySize(length)];
        }

        public BitArrayWorker(CompoundSet<C> compoundSet, int[] sequence) {
            this.compoundSet = compoundSet;
            this.sequence = sequence;
            this.length = sequence.length;
        }

        /**
         * This method should return the bit mask to be used to extract the
         * bytes you are interested in working with. See solid implementations
         * on how to create these
         */
        protected abstract byte bitMask();

        /**
         * Should return the maximum amount of compounds we can encode per int
         */
        protected abstract int compoundsPerDatatype();

        /**
         * Should return the inverse information that {@link #generateCompoundsToIndex() }
         * returns i.e. if the Compound C returns 1 from compoundsToIndex then we
         * should find that compound here in position 1
         */
        protected abstract List<C> generateIndexToCompounds();

        /**
         * Returns what the value of a compound is in the backing bit storage i.e.
         * in 2bit storage the value 0 is encoded as 00 (in binary).
         */
        protected abstract Map<C, Integer> generateCompoundsToIndex();

        /**
         * Returns how many bits are used to represent a compound e.g. 2 if using
         * 2bit encoding.
         */
        protected int bitsPerCompound() {
            return BYTES_PER_INT / compoundsPerDatatype();
        }

        public int seqArraySize(int length) {
            return (int) Math.ceil((double) length / (double) compoundsPerDatatype());
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
            hashcode = null;
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
            int masked = (int) (shifted & bitMask());

            //If we could encode 4 compounds then our max masked value is 3
            if (masked > (compoundsPerDatatype() - 1)) {
                throw new IllegalStateException("Got a masked value of " + masked + "; do not understand values greater than " + (compoundsPerDatatype() - 1));
            }
            return getIndexToCompoundsLookup().get(masked);
        }

        /**
         * Since bit encoding only supports a finite number of bases
         * it is more than likely when processing sequence you will encounter a
         * compound which is not covered by the encoding e.g. N in a 2bit sequence.
         * You can override this to convert the unknown base into one you can
         * process or store locations of unknown bases for a level of post processing
         * in your subclass.
         *
         * @param compound Compound process
         * @return Byte representation of the compound
         * @throws IllegalStateException Done whenever this method is invoked
         */
        protected byte processUnknownCompound(C compound, int position) throws IllegalStateException {
            throw new IllegalStateException("Do not know how to translate the compound " + compound + " to a " + bitsPerCompound() + "bit representation");
        }

        /**
         * Returns a list of compounds the index position of which is used
         * to translate from the byte representation into a compound.
         */
        protected List<C> getIndexToCompoundsLookup() {
            if (indexToCompoundsLookup == null) {
                indexToCompoundsLookup = generateIndexToCompounds();
            }
            return indexToCompoundsLookup;
        }

        /**
         * Returns a map which converts from compound to an integer representation
         */
        protected Map<C, Integer> getCompoundsToIndexLookup() {
            if (compoundsToIndexLookup == null) {
                compoundsToIndexLookup = generateCompoundsToIndex();
            }
            return compoundsToIndexLookup;
        }

        /**
         * Converting a biological index to the int which is used to store that
         * position's data.
         *
         * <ul>
         * <li>Assuming 2bit encoding using the 17 base in a sequence</li>
         * <li>Convert into 0 index; 17 - 1 = 16</li>
         * <li>Divide by the number of compounds per int; 16 / 16</li>
         * <li>Floor the value; floor(1) = 1</li>
         * </ul>
         *
         * Running this for position 13
         *
         * <ul>
         * <li>13 - 1 = 12</li>
         * <li>12 / 16 = 0.75</li>
         * <li>floor(0.75) = 0</li>
         * </ul>
         */
        private int biologicalIndexToArrayIndex(int index) {
            return ((index - 1) / compoundsPerDatatype());
        }

        /**
         * Convert from bio to 0 index, remainder & then multiply by 2
         * <ul>
         * <li>Using 2bit encoding and working with position 19</li>
         * <li>19 is the 3rd position in the second int</li>
         * <li>Means a shift of 4 into that int to get the right data out</li>
         * <li>Also must convert into the non-bio index</li>
         * <li>19 - 1 = 18</li>
         * <li>18 % compoundsPerDatatype() (16) = 2</li>
         * <li>2 * bits per compound (2) = 4</li>
         * </ul>
         */
        private byte shiftBy(int index) {
            return (byte) (((index - 1) % compoundsPerDatatype()) * bitsPerCompound());
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

        @Override
        public int hashCode() {
            if(hashcode == null) {
                int s = Hashcoder.SEED;
                s = Hashcoder.hash(s, sequence);
                s = Hashcoder.hash(s, indexToCompoundsLookup);
                s = Hashcoder.hash(s, compoundSet);
                hashcode = s;
            }
            return hashcode;
        }

        @Override
        @SuppressWarnings("unchecked")
        public boolean equals(Object o) {
            if(Equals.classEqual(this, o)) {
                BitArrayWorker<C> that = (BitArrayWorker<C>)o;
                return  Equals.equal(compoundSet, that.compoundSet) &&
                        Equals.equal(indexToCompoundsLookup, that.indexToCompoundsLookup) &&
                        Equals.equal(sequence, that.sequence);
            }
            return false;
        }
    }
}
