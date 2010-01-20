/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence;

import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.loader.SequenceStringProxyLoader;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.SequenceProxyLoader;

/**
 *
 * @author Scooter
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
