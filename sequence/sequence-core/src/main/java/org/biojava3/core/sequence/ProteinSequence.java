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

import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.loader.SequenceStringProxyLoader;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.SequenceProxyLoader;

/**
 *
 * @author Scooter Willis
 */
public class ProteinSequence extends AbstractSequence<AminoAcidCompound> {

    private DNASequence parentDNASequence = null;
    private Integer begin;
    private Integer end;

    public ProteinSequence(String seqString) {
        super(seqString, new AminoAcidCompoundSet());
    }

    public ProteinSequence(SequenceProxyLoader<AminoAcidCompound> proxyLoader) {
        super(proxyLoader, new AminoAcidCompoundSet());
    }

    public void setParentDNASequence(DNASequence parentDNASequence, Integer begin, Integer end) {
        this.begin = begin;
        this.end = end;
    }

    /**
     * @return the parentTranscriptSequence
     */
    public DNASequence getParentDNASequence() {
        return parentDNASequence;
    }

    /**
     * @return the begin
     */
    public Integer getBegin() {
        return begin;
    }

    /**
     * @return the end
     */
    public Integer getEnd() {
        return end;
    }

    public static void main(String[] args) {
        ProteinSequence proteinSequence = new ProteinSequence("ARNDCEQGHILKMFPSTWYVBZJX");
        System.out.println(proteinSequence.toString());

        SequenceStringProxyLoader<AminoAcidCompound> sequenceStringProxyLoader = new SequenceStringProxyLoader("XRNDCEQGHILKMFPSTWYVBZJA", AminoAcidCompoundSet.getAminoAcidCompoundSet());
        ProteinSequence proteinSequenceFromProxy = new ProteinSequence(sequenceStringProxyLoader);
        System.out.println(proteinSequenceFromProxy.toString());

    }
}
