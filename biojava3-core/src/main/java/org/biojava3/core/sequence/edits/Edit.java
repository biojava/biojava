package org.biojava3.core.sequence.edits;

import java.util.ArrayList;
import java.util.List;

import org.biojava3.core.sequence.BasicSequence;
import org.biojava3.core.sequence.storage.JoiningSequenceReader;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

/**
 * Interface for carrying out edit operations on a Sequence. The 3 major
 * methods of Editing are supported
 *
 * <ul>
 * <li>Insertion</li>
 * <li>Deletion</li>
 * <li>Substitution</li>
 * </ul>
 *
 * The interface is provided so end users can use our implementations, which
 * are implementations which attempts to create views of Sequences in an
 * editted form not a full-realised editted Sequence, or their own.
 *
 * @author ayates
 * @param <C> The type of compound to edit
 */
public interface Edit<C extends Compound> {

    Sequence<C> edit(Sequence<C> sequence);

    /**
     * Abstract class which defines all edit operations as a call to discover
     * what 5' and 3' ends of an editing Sequence should be joined together
     * with a target Sequence. These ends can be of 0 length but conceptionally
     * they can still exist.
     */
    public static abstract class AbstractEdit<C extends Compound> implements Edit<C> {

        /**
         * Should return the 5-prime end of the given Sequence according to
         * the edit. An empty Sequence is valid.
         */
        protected abstract Sequence<C> getFivePrime(Sequence<C> editingSequence);

        /**
         * Should return the 3-prime end of the given Sequence according to
         * the edit. An empty Sequence is valid.
         */
        protected abstract Sequence<C> getThreePrime(Sequence<C> editingSequence);

      
        public Sequence<C> edit(Sequence<C> editingSequence) {
            Sequence<C> targetSequence = getTargetSequence(editingSequence);
            List<Sequence<C>> sequences = new ArrayList<Sequence<C>>();

            sequences.add(getFivePrime(editingSequence));
            sequences.add(targetSequence);
            sequences.add(getThreePrime(editingSequence));

            return new JoiningSequenceReader<C>(sequences);
        }
        private int start = -1;
        private int end = -1;
        private String stringSequence;
        private Sequence<C> sequence;

        public AbstractEdit(int start) {
            this.start = start;
        }

        public AbstractEdit(int start, int end) {
            this.start = start;
            this.end = end;
        }

        protected void setStringSequence(String stringSequence) {
            this.stringSequence = stringSequence;
        }

        protected void setSequence(Sequence<C> sequence) {
            this.sequence = sequence;
        }

        /**
         * Returns the Sequence which is our edit.
         *
         * @param editingSequence Asked for in-case we need to do String to
         * Sequence conversion so we need a CompoundSet which is given
         * by the Sequence we are editing
         * @return The Sequence<C> object we wish to insert
         */
        public Sequence<C> getTargetSequence(Sequence<C> editingSequence) {
            if (sequence == null && stringSequence != null) {
                sequence = new BasicSequence<C>(
                        stringSequence, editingSequence.getCompoundSet());
            }
            return sequence;
        }

        /**
         * Returns an empty sequence with the given compound set of the editing
         * sequence
         */
        protected Sequence<C> getEmptySequence(Sequence<C> editingSequence) {
            return new BasicSequence<C>("", editingSequence.getCompoundSet());
        }

        public int getStart() {
            return start;
        }

        public int getEnd() {
            return end;
        }
    }

    /**
     * Implementation which allows for the deletion of bases from a Sequence
     */
    public static class Delete<C extends Compound> extends AbstractEdit<C> {

        public Delete(int position) {
            this(position, position);
        }

        public Delete(int start, int end) {
            super(start, end);
            setStringSequence("");
        }

        protected int getRealStart() {
            return getStart() - 1;
        }

        protected int getRealEnd() {
            return getEnd() + 1;
        }

        @Override
        protected Sequence<C> getFivePrime(Sequence<C> editingSequence) {
            int start = getRealStart();
            if (start == 0) {
                return getEmptySequence(editingSequence);
            }
            return editingSequence.getSubSequence(1, start);
        }

        @Override
        protected Sequence<C> getThreePrime(Sequence<C> editingSequence) {
            int end = getRealEnd();
            if (end > editingSequence.getLength()) {
                return getEmptySequence(editingSequence);
            }
            return editingSequence.getSubSequence(end, editingSequence.getLength());
        }
    }

    /**
     * Edit implementation which allows us to insert a base at any position
     * in a Sequence. Specifying 1 base is used to insert at the start and
     * end of a Sequence. If you wish to carry out an in-sequence insertion
     * then you specify the flanking base positions e.g.
     *
     * <pre>
     *   ACTG insert TT @ position 1   : TTACGT
     *   ACTG insert TT @ position 2,3 : ACTTGT
     *   ACTG insert A  @ position 4   : ACGTA
     * </pre>
     *
     * The code will raise exceptions if you attempt a single base edit
     * with an insertion.
     */
    public static class Insert<C extends Compound> extends AbstractEdit<C> {

        private final boolean singlePosition;

        public Insert(String sequence, int position) {
            super(position, position);
            this.singlePosition = true;
            setStringSequence(sequence);
        }

        public Insert(Sequence<C> sequence, int position) {
            super(position, position);
            this.singlePosition = true;
            setSequence(sequence);
        }

        public Insert(String sequence, int start, int stop) {
            super(start, stop);
            this.singlePosition = false;
            setStringSequence(sequence);
        }

        public Insert(Sequence<C> sequence, int start, int stop) {
            super(start, stop);
            this.singlePosition = false;
            setSequence(sequence);
        }

        @Override
        protected Sequence<C> getFivePrime(Sequence<C> editingSequence) {
            if (singlePosition) {
                if (getStart() == 1) {
                    return getEmptySequence(editingSequence);
                } else if (getEnd() == editingSequence.getLength()) {
                    return editingSequence;
                } else {
                    throw new IllegalStateException("Given one position to "
                            + "insert at but this is not the start or end "
                            + "of the Sequence; cannot support this");
                }
            }
            return editingSequence.getSubSequence(1, getStart());
        }

        @Override
        protected Sequence<C> getThreePrime(Sequence<C> editingSequence) {
            if (singlePosition) {
                if (getStart() == 1) {
                    return editingSequence;
                } else if (getEnd() == editingSequence.getLength()) {
                    return getEmptySequence(editingSequence);
                } else {
                    throw new IllegalStateException("Given one position to "
                            + "insert at but this is not the start or end "
                            + "of the Sequence; cannot support this");
                }
            }
            return editingSequence.getSubSequence(getEnd(), editingSequence.getLength());
        }
    }

    /**
     * Allows for the substitution of bases into an existing Sequence. This
     * allows us to do edits like:
     *
     * <pre>
     *    Sub TT @ position 2
     *    AAAA -> ATTA
     * </pre>
     *
     * We do not support
     *
     * Edits do not require the length of the insertion but do rely on the
     * presence of a CompoundSet to parse a String (if given) which means
     * the eventual length of a Sequence is a lazy operation.
     */
    public static class Substitute<C extends Compound> extends AbstractEdit<C> {

        public Substitute(String sequence, int position) {
            super(position);
            setStringSequence(sequence);
        }

        public Substitute(Sequence<C> sequence, int position) {
            super(position);
            setSequence(sequence);
        }

        /**
         * Must use this rather than the no-args getEnd as this can return
         * -1 and the length of a sub is dependent on the length of the
         * Sequence; we cannot assume 1:1 mapping between characters in a
         * String and the number of compounds we will have to insert.
         */
        public int getEnd(Sequence<C> sequence) {
            if (getEnd() == -1) {
                int start = getStart();
                int length = getTargetSequence(sequence).getLength();
                return (start + length) - 1;
            }
            return getEnd();
        }

        @Override
        protected Sequence<C> getFivePrime(Sequence<C> editingSequence) {
            int start = getStart();
            if (start == 1) {
                return getEmptySequence(editingSequence);
            }
            return editingSequence.getSubSequence(1, start - 1);
        }

        @Override
        protected Sequence<C> getThreePrime(Sequence<C> editingSequence) {
            int end = getEnd(editingSequence);
            if (end > editingSequence.getLength()) {
                throw new IndexOutOfBoundsException(end +
                        " is greater than the max index of " +
                        "the editing sequence (" +
                        editingSequence.getLength());
            } else if (end == editingSequence.getLength()) {
                return getEmptySequence(editingSequence);
            }
            return editingSequence.getSubSequence(end + 1, editingSequence.getLength());
        }
    }
}
