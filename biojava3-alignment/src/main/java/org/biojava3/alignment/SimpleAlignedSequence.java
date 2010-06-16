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

package org.biojava3.alignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import org.biojava3.alignment.template.AlignedSequence;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.location.SimpleLocation;
import org.biojava3.core.sequence.location.template.Location;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.template.SequenceView;

/**
 * Implements a data structure for a {@link Sequence} within an alignment.
 *
 * @author Mark Chapman
 * @param <C> each element of the {@link Sequence} is a {@link Compound} of type C
 */
public class SimpleAlignedSequence<C extends Compound> implements AlignedSequence<C> {

    private static final String gap = "-";

    // always stored
    private Sequence<C> original;
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
    public SimpleAlignedSequence(Sequence<C> original, List<Step> steps) {
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
    public SimpleAlignedSequence(Sequence<C> original, List<Step> steps, int numBefore, int numAfter) {
        this.original = original;
        this.numBefore = numBefore;
        this.numAfter = numAfter;
        length = steps.size();
        setLocation(steps);
    }

    // methods for AlignedSequence

    @Override
    public int getAlignmentIndexAt(int sequenceIndex) {
        if (alignmentFromSequence == null) {
            setAlignmentFromSequence();
        }
        return alignmentFromSequence[sequenceIndex - 1];
    }

    @Override
    public int getEnd() {
        return location.getEnd();
    }

    @Override
    public Location getLocationInAlignment() {
        return location;
    }

    @Override
    public int getNumGaps() {
        if (numGaps == -1) {
            setNumGaps();
        }
        return numGaps;
    }

    @Override
    public Sequence<C> getOriginalSequence() {
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
            setSequenceFromAlignment();
        }
        return sequenceFromAlignment[alignmentIndex - 1];
    }

    @Override
    public int getStart() {
        return location.getStart();
    }

    @Override
    public boolean isCircular() {
        return location.isCircular();
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
        return getSequenceAsString(1, length, Strand.UNDEFINED);
    }

    @Override
    public String getSequenceAsString(Integer start, Integer end, Strand strand) {
        StringBuilder s = new StringBuilder();
        for (C compound : getAsList().subList(start - 1, end)) {
            s.append(compound == null ? gap : getCompoundSet().getStringForCompound(compound));
        }
        return s.toString();
    }

    @Override
    public SequenceView<C> getSubSequence(Integer start, Integer end) {
        // TODO SequenceView<C> getSubSequence(Integer start, Integer end)
        return null;
    }

    // method for Iterable

    @Override
    public Iterator<C> iterator() {
        return getAsList().iterator();
    }

    /**
     * Provides standard Java language access to results of {@link #getSequenceAsString()}.
     */
    @Override
    public String toString() {
        return getSequenceAsString();
    }

    // helper methods

    // determines if this sequence has a gap at a particular alignment location
    private boolean isGap(int alignmentIndex) {
        for (Location sublocation : location) {
            if (sublocation.getStart() <= alignmentIndex && alignmentIndex <= sublocation.getEnd()) {
                return false;
            }
        }
        return true;
    }

    // initializes alignmentFromSequence
    private void setAlignmentFromSequence() {
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

    // initializes location
    private void setLocation(List<Step> steps) {
        List<Location> sublocations = new ArrayList<Location>();
        int step = 0, i = numBefore+numAfter, iMax = original.getLength();
        while (step < length && steps.get(step) == Step.GAP) {
            step++;
        }
        while (step < length && i < iMax) {
            int start = step + 1;
            while (step < length && steps.get(step) == Step.COMPOUND && i < iMax) {
                step++;
                i++;
            }
            sublocations.add(new SimpleLocation(start, step, Strand.UNDEFINED));
            while (step < length && steps.get(step) == Step.GAP) {
                step++;
            }
        }
        int start = sublocations.get(0).getStart(), end = sublocations.get(sublocations.size() - 1).getEnd();
        location = (sublocations.size() > 1) ? new SimpleLocation(start, end, Strand.UNDEFINED, sublocations) :
                new SimpleLocation(start, end, Strand.UNDEFINED);

        if (step != length || i != iMax) {
            throw new IllegalArgumentException("Given sequence does not fit in alignment.");
        }
    }

    // initializes the number of gaps
    private void setNumGaps() {
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

    // initializes sequenceFromAlignment
    private void setSequenceFromAlignment() {
        sequenceFromAlignment = new int[length];
        int a = 1, s = numBefore + 1;
        for (int i = 0; i < getStart(); i++, a++) {
            sequenceFromAlignment[a - 1] = s;
        }
        for (; a <= length; a++) {
            if (s < getEnd() && !isGap(a)) {
                s++;
            }
            sequenceFromAlignment[a - 1] = s;
        }
    }

}
