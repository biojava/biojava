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
import java.util.Collections;
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

    private Sequence<C> original;
    private List<Step> steps;
    private Location location;
    private C gap;
    private int gaps;
    private int[] alignmentFromSequence, sequenceFromAlignment;

    /**
     * Creates an {@link AlignedSequence} for the given {@link Sequence}.
     *
     * @param original the original {@link Sequence} before alignment
     * @param steps lists whether the sequence aligns a {@link Compound} or gap at each index of the alignment
     */
    public SimpleAlignedSequence(Sequence<C> original, List<Step> steps) {
        this.original = original;
        this.steps = Collections.unmodifiableList(steps);
        List<Location> sublocations = new ArrayList<Location>();
        int step = 0;
        while (step < steps.size() && steps.get(step) == Step.GAP) {
            step++;
        }
        for (int i = step, start; i < original.getLength(); i++) {
            start = step;
            while (step < steps.size() && steps.get(step) == Step.COMPOUND) {
                step++;
                i++;
            }
            sublocations.add(new SimpleLocation(start + 1, step, Strand.UNDEFINED));
            while (step < steps.size() && steps.get(step) == Step.GAP) {
                step++;
            }
            if (step < steps.size()) {
                gaps++;
            }
        }
        int start = sublocations.get(0).getStart(), end = sublocations.get(sublocations.size() - 1).getEnd();
        location = (sublocations.size() > 1) ? new SimpleLocation(start, end, Strand.UNDEFINED, sublocations) :
                new SimpleLocation(start, end, Strand.UNDEFINED);
        gap = original.getCompoundSet().getCompoundForString("-");
    }

    @Override
    public int getAlignmentIndexAt(int sequenceIndex) {
        if (alignmentFromSequence == null) {
            alignmentFromSequence = new int[original.getLength()];
            for (int i = 0, step = 0; i < alignmentFromSequence.length; i++) {
                while (step < steps.size() && steps.get(step) == Step.GAP) {
                    step++;
                }
                if (step < steps.size()) {
                    step++;
                }
                alignmentFromSequence[i] = step;
            }
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
        return gaps;
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
            sequenceFromAlignment = new int[steps.size()];
            for (int i = 0; i < sequenceFromAlignment.length; i++) {
                sequenceFromAlignment[i] = (i < location.getStart()) ? 1 : ((steps.get(i) == Step.GAP) ?
                        sequenceFromAlignment[i - 1] : sequenceFromAlignment[i - 1] + 1);
            }
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

    @Override
    public int countCompounds(C... compounds) {
        return original.countCompounds(compounds);
    }

    @Override
    public AccessionID getAccession() {
        return original.getAccession();
    }

    @Override
    public List<C> getAsList() {
        return original.getAsList();
    }

    @Override
    public C getCompoundAt(int position) {
        return (steps.get(position - 1) == Step.GAP) ? gap : original.getCompoundAt(getSequenceIndexAt(position));
    }

    @Override
    public CompoundSet<C> getCompoundSet() {
        return original.getCompoundSet();
    }

    @Override
    public int getIndexOf(C compound) {
        return getAlignmentIndexAt(original.getIndexOf(compound));
    }

    @Override
    public int getLastIndexOf(C compound) {
        return getAlignmentIndexAt(original.getLastIndexOf(compound));
    }

    @Override
    public int getLength() {
        return steps.size();
    }

    @Override
    public String getSequenceAsString() {
        return getSequenceAsString(1, steps.size(), Strand.UNDEFINED);
    }

    @Override
    public String getSequenceAsString(Integer start, Integer end, Strand strand) {
        StringBuilder s = new StringBuilder();
        for (int i = start - 1; i < end; i++) {
            C compound = getCompoundAt(i);
            s.append((compound == null) ? "-" : original.getCompoundSet().getStringForCompound(compound));
        }
        return s.toString();
    }

    @Override
    public SequenceView<C> getSubSequence(Integer start, Integer end) {
        // TODO SequenceView<C> getSubSequence(Integer start, Integer end)
        return null;
    }

    @Override
    public Iterator<C> iterator() {
        List<C> compounds = original.getAsList();
        for (int i = 0; i < steps.size(); i++) {
            if (steps.get(i) == Step.GAP) {
                compounds.add(i, gap);
            }
        }
        return compounds.iterator();
    }

}
