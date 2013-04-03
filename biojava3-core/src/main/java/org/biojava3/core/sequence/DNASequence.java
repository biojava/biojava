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
 * Created on DATE
 *
 */
package org.biojava3.core.sequence;

import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.loader.StringProxySequenceReader;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.SequenceMixin;
import org.biojava3.core.sequence.template.SequenceReader;
import org.biojava3.core.sequence.template.SequenceView;
import org.biojava3.core.sequence.transcription.Frame;
import org.biojava3.core.sequence.transcription.TranscriptionEngine;
import org.biojava3.core.sequence.views.ComplementSequenceView;
import org.biojava3.core.sequence.views.ReversedSequenceView;

/**
 * This is class should model the attributes associated with a DNA sequence
 *
 * @author Scooter Willis
 */
public class DNASequence extends AbstractSequence<NucleotideCompound> {
/**
 * The type of DNA sequence
 */
    public enum DNAType {
        CHROMOSOME, MITOCHONDRIAL, PLASMID, PLASTID, UNKNOWN
    }
    private DNAType dnaType = DNAType.UNKNOWN;

    /**
     * Shouldn't be used but makes it bean happy
     */
    public DNASequence() {
//        throw new UnsupportedOperationException("Null constructor not supported");
    }

    /**
     * String is king and create a sequence from DNA with default DNA compound set
     * @param seqString
     */
    public DNASequence(String seqString) {
        super(seqString, DNACompoundSet.getDNACompoundSet());
    }

    /**
     * Create a sequence where the actual storage of the sequence data is somewhere else
     * @param proxyLoader
     */
    public DNASequence(SequenceReader<NucleotideCompound> proxyLoader) {
        super(proxyLoader, DNACompoundSet.getDNACompoundSet());
    }

    /**
     * Create a sequence from a string with user defined compound set
     * @param seqString
     * @param compoundSet
     */
    public DNASequence(String seqString, CompoundSet<NucleotideCompound> compoundSet) {
        super(seqString, compoundSet);
    }

    /**
     * Create a sequence from a ProxySequencereader and user defined compound set
     * @param proxyLoader
     * @param compoundSet
     */
    public DNASequence(SequenceReader<NucleotideCompound> proxyLoader, CompoundSet<NucleotideCompound> compoundSet) {
        super(proxyLoader, compoundSet);
    }

    /**
     * Return the RNASequence equivalent of the DNASequence using default Transcription Engine. Not all
     * species follow the same rules. If you don't know better use this method
     * @return RNA sequence
     */
    public RNASequence getRNASequence() {
      return getRNASequence(Frame.getDefaultFrame());
    }

    /**
     * Allow a user to pass in a rules engine to do the DNA to RNA translation
     * @param engine
     * @return RNA sequence
     */ 
    public RNASequence getRNASequence(TranscriptionEngine engine) {
      return getRNASequence(engine, Frame.getDefaultFrame());
    }

    /**
     * Allows the user to pass in the Frame shift.
     * @param frame
     * @return rna sequence
     */
    public RNASequence getRNASequence(Frame frame) {
      return getRNASequence(TranscriptionEngine.getDefault(), frame);
    }

    public RNASequence getRNASequence(TranscriptionEngine engine, Frame frame) {
      return (RNASequence) engine.getDnaRnaTranslator().createSequence(this, frame);
    }

    /**
     * Get the GC count in the DNA Sequence
     * @return GC count
     */
    public int getGCCount() {
        return SequenceMixin.countGC(this);
    }

    /**
     * Returns a Sequence which runs in the current reverse order
     */
    public SequenceView<NucleotideCompound> getReverse() {
        return new ReversedSequenceView<NucleotideCompound>(this);
    }

    /**
     * Returns a Sequence which will complement every base
     */
    public SequenceView<NucleotideCompound> getComplement() {
        return new ComplementSequenceView<NucleotideCompound>(this);
    }

    /**
     * Delegates to {@link #getInverse() } for the reverse complement
     */
    public SequenceView<NucleotideCompound> getReverseComplement() {
        return getInverse();
    }

    /**
     * @return the dnaType
     */
    public DNAType getDNAType() {
        return dnaType;
    }

    /**
     * @param dnaType the dnaType to set
     */
    public void setDNAType(DNAType dnaType) {
        this.dnaType = dnaType;
    }

    public static void main(String[] args) {
        DNASequence dnaSequence = new DNASequence("ATCG");
        System.out.println(dnaSequence.toString());

        StringProxySequenceReader<NucleotideCompound> sequenceStringProxyLoader =
                new StringProxySequenceReader<NucleotideCompound>("GCTA", DNACompoundSet.getDNACompoundSet());
        DNASequence dnaSequenceFromProxy = new DNASequence(sequenceStringProxyLoader);
        System.out.println(dnaSequenceFromProxy.toString());

    }
}
