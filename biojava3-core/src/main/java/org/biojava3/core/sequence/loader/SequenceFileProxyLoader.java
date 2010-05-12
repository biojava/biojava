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

import java.io.File;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.biojava3.core.sequence.template.SequenceProxyView;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.exceptions.CompoundNotFoundError;
import org.biojava3.core.exceptions.FileAccessError;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.Strand;

import org.biojava3.core.sequence.io.template.SequenceParserInterface;
import org.biojava3.core.sequence.storage.SequenceAsStringHelper;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.ProxySequenceReader;
import org.biojava3.core.sequence.template.SequenceMixin;
import org.biojava3.core.sequence.template.SequenceView;

public class SequenceFileProxyLoader<C extends Compound> implements ProxySequenceReader<C> {

    SequenceParserInterface sequenceParser;
    private CompoundSet<C> compoundSet;
    private List<C> parsedCompounds = new ArrayList<C>();
    File file;
    long sequenceStartIndex = -1;
    int sequenceLength = -1;
    private boolean initialized = false;

    public SequenceFileProxyLoader(File file, SequenceParserInterface sequenceParser, long sequenceStartIndex, int sequenceLength, CompoundSet<C> compoundSet) {
        this.sequenceParser = sequenceParser;
        this.file = file;
        this.sequenceStartIndex = sequenceStartIndex;
        this.sequenceLength = sequenceLength;
        setCompoundSet(compoundSet);
    }

    public void setCompoundSet(CompoundSet<C> compoundSet) {
        this.compoundSet = compoundSet;
    }

    private boolean init() {
        try {
            RandomAccessFile randomAccessFile = new RandomAccessFile(file, "r");
            randomAccessFile.seek(sequenceStartIndex);
            String sequence = sequenceParser.getSequence(randomAccessFile, sequenceLength);
            setContents(sequence);
        } catch (Exception e) {
            throw new FileAccessError("Error accessing " + file + " offset=" + sequenceStartIndex + " sequenceLength=" + sequenceLength + " " + e.toString());
        }
        return true;
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

        setInitialized(true);
    }

    public int getLength() {
        return sequenceLength;
    }

    public C getCompoundAt(int position) {
        if (this.isInitialized()) {
            init();
        }
        return this.parsedCompounds.get(position - 1);
    }

    public int getIndexOf(C compound) {
        if (this.isInitialized()) {
            init();
        }
        return this.parsedCompounds.indexOf(compound) + 1;
    }

    public int getLastIndexOf(C compound) {
        if (this.isInitialized()) {
            init();
        }
        return this.parsedCompounds.lastIndexOf(compound) + 1;
    }

    @Override
    public String toString() {
        if (this.isInitialized()) {
            init();
        }
        return getSequenceAsString();
    }

    @Override
    public String getSequenceAsString() {
        return getSequenceAsString(1, getLength(),Strand.POSITIVE);
    }

    public String getSequenceAsString(Integer bioBegin, Integer bioEnd,Strand strand) {
        if (this.isInitialized()) {
            init();
        }
        SequenceAsStringHelper<C> sequenceAsStringHelper = new SequenceAsStringHelper<C>();
        return sequenceAsStringHelper.getSequenceAsString(this.parsedCompounds, compoundSet, bioBegin, bioEnd, strand);
    }

    public List<C> getAsList() {
        if (this.isInitialized()) {
            init();
        }
        return this.parsedCompounds;

    }

    public SequenceView<C> getSubSequence(final Integer bioBegin, final Integer bioEnd) {
        if (this.isInitialized()) {
            init();
        }
        return new SequenceProxyView<C>(SequenceFileProxyLoader.this,bioBegin, bioEnd);
    }

    public Iterator<C> iterator() {
        if (this.isInitialized()) {
            init();
        }
        return this.parsedCompounds.iterator();
    }

    public CompoundSet<C> getCompoundSet() {
        return compoundSet;
    }

    /**
     * @return the initialized
     */
    public boolean isInitialized() {
        return initialized;
    }

    /**
     * @param initialized the initialized to set
     */
    public void setInitialized(boolean initialized) {
        this.initialized = initialized;
    }

    @Override
    public AccessionID getAccession() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public int countCompounds(C... compounds) {
      return SequenceMixin.countCompounds(this, compounds);
    }
}
