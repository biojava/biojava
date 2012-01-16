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

/**
 * This class represents the storage container of a sequence stored in a fasta file where
 * the initial parsing of the file we store the offset and length of the sequence. When a call
 * is made to any method that needs sequence data then the file will be opened and the sequence
 * loaded. This class could be improved by using the hints or a some algorithm that indicates
 * the sequence data once loaded should stay loaded. Could keep track of the last time sequence
 * data was loaded and then after X amount of time clear the contents to free up memory.
 *
 *
 * @author Scooter Willis <willishf at gmail dot com>
 * @param <C>
 */
public class SequenceFileProxyLoader<C extends Compound> implements ProxySequenceReader<C> {

    SequenceParserInterface sequenceParser;
    private CompoundSet<C> compoundSet;
    private List<C> parsedCompounds = new ArrayList<C>();
    File file;
    long sequenceStartIndex = -1;
    int sequenceLength = -1;
    private boolean initialized = false;

    /**
     *
     * @param file The file where the sequence will be found
     * @param sequenceParser The parser to use to load the sequence
     * @param sequenceStartIndex The file offset to the start of the sequence
     * @param sequenceLength The length of the sequence
     * @param compoundSet
     */
    public SequenceFileProxyLoader(File file, SequenceParserInterface sequenceParser, long sequenceStartIndex, int sequenceLength, CompoundSet<C> compoundSet) {
        this.sequenceParser = sequenceParser;
        this.file = file;
        this.sequenceStartIndex = sequenceStartIndex;
        this.sequenceLength = sequenceLength;
        setCompoundSet(compoundSet);
    }

    /**
     *
     * @param compoundSet
     */
    public void setCompoundSet(CompoundSet<C> compoundSet) {
        this.compoundSet = compoundSet;
    }

    /**
     *  Load the sequence
     * @return
     */
    private boolean init() {
        try {
            RandomAccessFile randomAccessFile = new RandomAccessFile(file, "r");
            randomAccessFile.seek(sequenceStartIndex);
            String sequence = sequenceParser.getSequence(randomAccessFile, sequenceLength);
            setContents(sequence);
            randomAccessFile.close(); // close file to prevent too many being open
        } catch (Exception e) {
            throw new FileAccessError("Error accessing " + file + " offset=" + sequenceStartIndex + " sequenceLength=" + sequenceLength + " " + e.toString());
        }
        return true;
    }

    /**
     *
     * @param sequence
     */
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

    /**
     *
     * @return
     */
    public int getLength() {
        return sequenceLength;
    }

    /**
     *
     * @param position
     * @return
     */
    public C getCompoundAt(int position) {
        if (!this.isInitialized()) {
            init();
        }
        return this.parsedCompounds.get(position - 1);
    }

    /**
     *
     * @param compound
     * @return
     */
    public int getIndexOf(C compound) {
        if (!this.isInitialized()) {
            init();
        }
        return this.parsedCompounds.indexOf(compound) + 1;
    }

    /**
     *
     * @param compound
     * @return
     */
    public int getLastIndexOf(C compound) {
        if (!this.isInitialized()) {
            init();
        }
        return this.parsedCompounds.lastIndexOf(compound) + 1;
    }

    /**
     *
     * @return
     */
    public String toString() {
        if (!this.isInitialized()) {
            init();
        }
        return getSequenceAsString();
    }

    /**
     *
     * @return
     */
    public String getSequenceAsString() {
        return getSequenceAsString(1, getLength(), Strand.POSITIVE);
    }

    /**
     *
     * @param bioBegin
     * @param bioEnd
     * @param strand
     * @return
     */
    public String getSequenceAsString(Integer bioBegin, Integer bioEnd, Strand strand) {

        if (!this.isInitialized()) {
            init();
        }
        SequenceAsStringHelper<C> sequenceAsStringHelper = new SequenceAsStringHelper<C>();
        return sequenceAsStringHelper.getSequenceAsString(this.parsedCompounds, compoundSet, bioBegin, bioEnd, strand);
    }

    /**
     *
     * @return
     */
    public List<C> getAsList() {
        if (!this.isInitialized()) {
            init();
        }
        return this.parsedCompounds;

    }

    /**
     *
     * @param bioBegin
     * @param bioEnd
     * @return
     */
    public SequenceView<C> getSubSequence(final Integer bioBegin, final Integer bioEnd) {
        if (!this.isInitialized()) {
            init();
        }
        return new SequenceProxyView<C>(SequenceFileProxyLoader.this, bioBegin, bioEnd);
    }

    /**
     *
     * @return
     */
    public Iterator<C> iterator() {
        if (!this.isInitialized()) {
            init();
        }
        return this.parsedCompounds.iterator();
    }

    /**
     *
     * @return
     */
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

    /**
     *
     * @return
     */
    public AccessionID getAccession() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     *
     * @param compounds
     * @return
     */
    public int countCompounds(C... compounds) {
        return SequenceMixin.countCompounds(this, compounds);
    }

    /**
     *
     * @return
     */
    @Override
    public SequenceView<C> getInverse() {
        return SequenceMixin.inverse(this);
    }
}
