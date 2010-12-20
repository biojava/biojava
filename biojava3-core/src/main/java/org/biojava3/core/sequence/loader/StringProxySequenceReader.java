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
 * @auther Scooter Willis
 *
 */
package org.biojava3.core.sequence.loader;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import org.biojava3.core.sequence.AccessionID;

import org.biojava3.core.sequence.template.SequenceProxyView;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.exceptions.CompoundNotFoundError;
import org.biojava3.core.sequence.Strand;

import org.biojava3.core.sequence.storage.SequenceAsStringHelper;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.ProxySequenceReader;
import org.biojava3.core.sequence.template.SequenceMixin;
import org.biojava3.core.sequence.template.SequenceView;

/**
 * An example of a ProxySequenceReader that is created from a String. Used for testing
 * @author Scooter Willis <willishf at gmail dot com>
 * @param <C>
 */

public class StringProxySequenceReader<C extends Compound> implements ProxySequenceReader<C> {

    private String sequence;
    private CompoundSet<C> compoundSet;
    private List<C> parsedCompounds = new ArrayList<C>();
    

    public StringProxySequenceReader(String sequence, CompoundSet<C> compoundSet) {
        this.sequence = sequence;
        setCompoundSet(compoundSet);
        setContents(sequence);
    }

    public void setCompoundSet(CompoundSet<C> compoundSet) {
        this.compoundSet = compoundSet;
    }

    public void setContents(String sequence) {
        // Horrendously inefficient - pretty much the way the old BJ did things.
        // TODO Should be optimised.
        this.parsedCompounds.clear();
        for (int i = 0; i < sequence.length();) {
            String compoundStr = null;
            C compound = null;
            for (int compoundStrLength = 1; compound == null && compoundStrLength <= compoundSet.getMaxSingleCompoundStringLength(); compoundStrLength++) {
                compoundStr = sequence.substring(i, i + compoundStrLength);
                compound = compoundSet.getCompoundForString(compoundStr);
            }
            if (compound == null) {
                throw new CompoundNotFoundError(compoundStr);
            } else {
                i += compoundStr.length();
            }
            this.parsedCompounds.add(compound);
        }
    }

    public int getLength() {
        return this.parsedCompounds.size();
    }

    public C getCompoundAt(int position) {
        return this.parsedCompounds.get(position - 1);
    }

    public int getIndexOf(C compound) {
        return this.parsedCompounds.indexOf(compound) + 1;
    }

    public int getLastIndexOf(C compound) {
        return this.parsedCompounds.lastIndexOf(compound) + 1;
    }

    
    public String toString() {
        return getSequenceAsString();
    }

    public String getSequenceAsString() {
        return sequence;
    }

    public List<C> getAsList() {
        return this.parsedCompounds;
    }


        
    public String getSequenceAsString(Integer bioBegin, Integer bioEnd,Strand strand) {
        SequenceAsStringHelper<C> sequenceAsStringHelper = new SequenceAsStringHelper<C>();
        return sequenceAsStringHelper.getSequenceAsString(this.parsedCompounds, compoundSet, bioBegin, bioEnd, strand);
    }

    public SequenceView<C> getSubSequence(final Integer bioBegin, final Integer bioEnd) {
        return new SequenceProxyView<C>(StringProxySequenceReader.this,bioBegin,bioEnd);
    }

    public Iterator<C> iterator() {
        return this.parsedCompounds.iterator();
    }

    public CompoundSet<C> getCompoundSet() {
      return compoundSet;
    }

    
    public AccessionID getAccession() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    
    public int countCompounds(C... compounds) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public SequenceView<C> getInverse() {
        return SequenceMixin.inverse(this);
    }
}
