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
 * Created on 01-21-2010
 *
 * @author Richard Holland
 *
 *
 */
package org.biojava3.core.sequence.template;

import java.util.Iterator;
import java.util.List;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.Strand;

public class SequenceProxyView<C extends Compound> implements SequenceView<C> {

    private Integer bioStart;
    private Integer bioEnd;
    private Sequence<C> sequence;

    public SequenceProxyView() {
    }

    public SequenceProxyView(Sequence<C> sequence) {
        this(sequence, 1, sequence.getLength());
    }

    public SequenceProxyView(Sequence<C> sequence, Integer bioStart, Integer bioEnd) {
        this.sequence = sequence;
        this.bioStart = bioStart;
        this.bioEnd = bioEnd;
    }

  
    public Sequence<C> getViewedSequence() {
        return sequence;
    }


    public String getSequenceAsString(Integer start, Integer end, Strand strand) {
        return getViewedSequence().getSequenceAsString(start, end, strand);
    }


    public String getSequenceAsString() {
        // TODO Optimise.
        return getSequenceAsString(getBioStart(), getBioEnd(), Strand.POSITIVE);
    }

    
    public List<C> getAsList() {
        return SequenceMixin.toList(this);
    }

    
    public C getCompoundAt(int position) {
        return getViewedSequence().getCompoundAt((getBioStart() + position) - 1);
    }

    
    public int getIndexOf(C compound) {
        return (getViewedSequence().getIndexOf(compound) + getBioStart()) - 1;
    }

    
    public int getLastIndexOf(C compound) {
        return (getViewedSequence().getLastIndexOf(compound) + getBioStart()) - 1;
    }

    
    public int getLength() {
        return (getBioEnd() - getBioStart()) + 1;
    }

    
    public CompoundSet<C> getCompoundSet() {
        return getViewedSequence().getCompoundSet();
    }

    
    public SequenceView<C> getSubSequence(final Integer bioStart, final Integer bioEnd) {
        return new SequenceProxyView<C>(this, bioStart, bioEnd);
    }

    
    public Iterator<C> iterator() {
        return new SequenceMixin.SequenceIterator<C>(this);
    }

    
    public AccessionID getAccession() {
        return getViewedSequence().getAccession();
    }

    /**
     * @return the bioStart
     */
    
    public Integer getBioStart() {
        return bioStart;
    }

    /**
     * @param bioStart the bioStart to set
     */
    public void setBioStart(Integer bioStart) {
        this.bioStart = bioStart;
    }

    /**
     * @return the bioEnd
     */
    
    public Integer getBioEnd() {
        return bioEnd;
    }

    /**
     * @param bioEnd the bioEnd to set
     */
    public void setBioEnd(Integer bioEnd) {
        this.bioEnd = bioEnd;
    }

    
    public int countCompounds(C... compounds) {
        return SequenceMixin.countCompounds(this, compounds);
    }
}
