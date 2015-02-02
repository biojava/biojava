/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on June 14, 2010
 * Author: Mark Chapman
 */

package org.biojava.nbio.alignment;

import org.biojava.nbio.alignment.template.AlignedSequence;
import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.Strand;
import org.biojava.nbio.core.sequence.location.SimpleLocation;
import org.biojava.nbio.core.sequence.location.template.Location;
import org.biojava.nbio.core.sequence.location.template.Point;
import org.biojava.nbio.core.sequence.template.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/**
 * Implements a data structure for a {@link Sequence} within an alignment.
 *
 * @author Mark Chapman
 * @param <C> each element of the {@link Sequence} is a {@link Compound} of type C
 */
public class SimpleAlignedSequence<S extends Sequence<C>, C extends Compound> implements AlignedSequence<S, C> {

    private static final String gap = "-";

    // always stored
    private AlignedSequence<S, C> prev;
    private S original;
    private int length, numBefore, numAfter;
    private Location location;

    // cached (lazily initialized)
    private int numGaps = -1;
    private int[] alignmentFromSequence, sequenceFromAlignment;

    /**
     * Creates an {@link AlignedSequence} for the given {@link Sequence} in a global alignment.
     *
     * @param original the original {@link Sequence} before alignment
     * @param steps lists whether the sequence aligns a {@link Compound} or gap at each index of the alignment
     * @throws IllegalArgumentException if given sequence does not fit in alignment
     */
    public SimpleAlignedSequence(S original, List<Step> steps) {
        this(original, steps, 0, 0);
    }

    /**
     * Creates an {@link AlignedSequence} for the given {@link Sequence} in a local alignment.
     *
     * @param original the original {@link Sequence} before alignment
     * @param steps lists whether the sequence aligns a {@link Compound} or gap at each index of the alignment
     * @param numBefore number of {@link Compound}s before a local alignment
     * @param numAfter number of {@link Compound}s after a local alignment
     * @throws IllegalArgumentException if given sequence does not fit in alignment
     */
    public SimpleAlignedSequence(S original, List<Step> steps, int numBefore, int numAfter) {
        this.original = original;
        this.numBefore = numBefore;
        this.numAfter = numAfter;
        length = steps.size();
        setLocation(steps);
    }

    /**
     * Creates a new {@link AlignedSequence} for the given {@link AlignedSequence} in a global alignment.
     *
     * @param prev the previous {@link AlignedSequence} before this alignment
     * @param steps lists whether the sequence aligns a {@link Compound} or gap at each index of the alignment
     * @throws IllegalArgumentException if given sequence does not fit in alignment
     */
    public SimpleAlignedSequence(AlignedSequence<S, C> prev, List<Step> steps) {
        this(prev, steps, 0, 0);
    }

    /**
     * Creates a new {@link AlignedSequence} for the given {@link AlignedSequence} in a local alignment.
     *
     * @param prev the previous {@link AlignedSequence} before this alignment
     * @param steps lists whether the sequence aligns a {@link Compound} or gap at each index of the alignment
     * @param numBefore number of {@link Compound}s before a local alignment
     * @param numAfter number of {@link Compound}s after a local alignment
     * @throws IllegalArgumentException if given sequence does not fit in alignment
     */
    public SimpleAlignedSequence(AlignedSequence<S, C> prev, List<Step> steps, int numBefore, int numAfter) {
        this.prev = prev;
        this.original = prev.getOriginalSequence();
        this.numBefore = numBefore;
        this.numAfter = numAfter;
        if (prev instanceof SimpleAlignedSequence<?, ?>) {
            SimpleAlignedSequence<?, ?> p = (SimpleAlignedSequence<?, ?>) prev;
            this.numBefore += p.numBefore;
            this.numAfter += p.numAfter;
        }
        length = steps.size();
        setLocation(steps);
    }

    // methods for AlignedSequence

    @Override
    public void clearCache() {
        alignmentFromSequence = null;
        sequenceFromAlignment = null;
    }

    @Override
    public int getAlignmentIndexAt(int sequenceIndex) {
        if (alignmentFromSequence == null) {
            alignmentFromSequence = new int[original.getLength()];
            int s = 1, a = 1;
            for (int i = 0; i < numBefore; i++, s++) {
                alignmentFromSequence[s - 1] = a;
            }
            for (; s <= alignmentFromSequence.length && a <= length; s++, a++) {
                while (a <= length && isGap(a)) {
                    a++;
                }
                alignmentFromSequence[s - 1] = a;
            }
            a--;
            for (int i = 0; i < numAfter; i++, s++) {
                alignmentFromSequence[s - 1] = a;
            }
        }
        return alignmentFromSequence[sequenceIndex - 1];
    }

    @Override
    public Point getEnd() {
        return location.getEnd();
    }

    @Override
    public Location getLocationInAlignment() {
        return location;
    }

    @Override
    public int getNumGaps() {
        if (numGaps == -1) {
            numGaps = 0;
            C cGap = getCompoundSet().getCompoundForString(gap);
            boolean inGap = false;
            for (C compound : getAsList()) {
                if (compound == null || compound.equalsIgnoreCase(cGap)) {
                    if (!inGap) {
                        numGaps++;
                        inGap = true;
                    }
                } else {
                    inGap = false;
                }
            }
        }
        return numGaps;
    }

    @Override
    public S getOriginalSequence() {
        return original;
    }

    @Override
    public int getOverlapCount() {
        // TODO handle circular alignments
        return 1;
    }

    @Override
    public int getSequenceIndexAt(int alignmentIndex) {
        if (sequenceFromAlignment == null) {
            sequenceFromAlignment = new int[length];
            int a = 1, s = numBefore + 1;
            for (int i = 0; i < getStart().getPosition(); i++, a++) {
                sequenceFromAlignment[a - 1] = s;
            }
            for (; a <= length; a++) {
                if (!isGap(a)) {
                    s++;
                }
                sequenceFromAlignment[a - 1] = s;
            }
        }
        return sequenceFromAlignment[alignmentIndex - 1];
    }

    @Override
    public Point getStart() {
        return location.getStart();
    }

    @Override
    public boolean isCircular() {
        return location.isCircular();
    }

    @Override
    public boolean isGap(int alignmentIndex) {
        if (getStart().getPosition() <= alignmentIndex && alignmentIndex <= getEnd().getPosition()) {
            if (!location.isComplex()) {
                return false;
            }
            for (Location sublocation : location) {
                if (sublocation.getStart().getPosition() <= alignmentIndex &&
                        alignmentIndex <= sublocation.getEnd().getPosition()) {
                    return false;
                }
            }
        }
        return true;
    }

    // methods for Sequence

    @Override
    public int countCompounds(C... compounds) {
        int count = 0;
        List<C> search = Arrays.asList(compounds);
        for (C compound : getAsList()) {
            if (search.contains(compound)) {
                count++;
            }
        }
        return count;
    }

    @Override
    public AccessionID getAccession() {
        return original.getAccession();
    }

    @Override
    public List<C> getAsList() {
        List<C> compounds = new ArrayList<C>();
        for (int i = 1; i <= length; i++) {
            compounds.add(getCompoundAt(i));
        }
        return compounds;
    }

    @Override
    public C getCompoundAt(int alignmentIndex) {
        return alignmentIndex >= 1 && alignmentIndex <= length && isGap(alignmentIndex) ?
                getCompoundSet().getCompoundForString(gap) :
                original.getCompoundAt(getSequenceIndexAt(alignmentIndex));
    }

    @Override
    public CompoundSet<C> getCompoundSet() {
        return original.getCompoundSet();
    }

    @Override
    public int getIndexOf(C compound) {
        for (int i = 1; i <= length; i++) {
            if (compound.equals(getCompoundAt(i))) {
                return i;
            }
        }
        return -1;
    }

    @Override
    public int getLastIndexOf(C compound) {
        for (int i = length; i >= 1; i--) {
            if (compound.equals(getCompoundAt(i))) {
                return i;
            }
        }
        return -1;
    }

    @Override
    public int getLength() {
        return length;
    }

    @Override
    public String getSequenceAsString() {
        return SequenceMixin.toString(this);
    }

    @Override
    public SequenceView<C> getSubSequence(Integer start, Integer end) {
        return SequenceMixin.createSubSequence(this, start, end);
    }

    // method for Iterable

    @Override
    public Iterator<C> iterator() {
        return getAsList().iterator();
    }

    // method from Object

    /**
     * Provides standard Java language access to results of {@link #getSequenceAsString()}.
     */
    @Override
    public String toString() {
        return getSequenceAsString();
    }

    // helper method to initialize the location
    private void setLocation(List<Step> steps) {
        List<Location> sublocations = new ArrayList<Location>();
        int start = 0, step = 0, oStep = numBefore+numAfter, oMax = this.original.getLength(), pStep = 0, pMax =
                (prev == null) ? 0 : prev.getLength();
        boolean inGap = true;

        // build sublocations: pieces of sequence separated by gaps
        for (; step < length; step++) {
            boolean isGapStep = (steps.get(step) == Step.GAP),
                    isGapPrev = (pStep < pMax && prev.isGap(pStep + 1));
            if (!isGapStep && !isGapPrev) {
                oStep++;
                if (inGap) {
                    inGap = false;
                    start = step + 1;
                }
            } else if (!inGap) {
                inGap = true;
                sublocations.add(new SimpleLocation(start, step, Strand.UNDEFINED));
            }
            if (prev != null && !isGapStep) {
                pStep++;
            }
        }
        if (!inGap) {
            sublocations.add(new SimpleLocation(start, step, Strand.UNDEFINED));
        }

        // combine sublocations into 1 Location
        if (sublocations.size() == 0) {
        	location = null;
        } else if (sublocations.size() == 1) {
        	location = sublocations.get(0);
        } else {
        	location = new SimpleLocation(sublocations.get(0).getStart(), sublocations.get(sublocations.size() - 1).getEnd(),
        			Strand.UNDEFINED,
                    false, sublocations);
        }
        // TODO handle circular alignments

        // check that alignment has correct number of compounds in it to fit original sequence
        if (step != length || oStep != oMax || pStep != pMax) {
            throw new IllegalArgumentException("Given sequence does not fit in alignment.");
        }
    }

    @Override
    //TODO Needs to implements
    public SequenceView<C> getInverse() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
