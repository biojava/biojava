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
import org.biojava3.core.sequence.loader.StringProxySequenceReader;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.ProxySequenceReader;

/**
 * The representation of a ProteinSequence
 * @author Scooter Willis
 */
public class ProteinSequence extends AbstractSequence<AminoAcidCompound> {

    /**
     * Create a protein from a string
     * @param seqString
     */
    public ProteinSequence(String seqString) {
        this(seqString, AminoAcidCompoundSet.getAminoAcidCompoundSet());
    }

    /**
     * Create a protein from a string with a user defined set of amino acids
     * @param seqString
     * @param compoundSet
     */
    public ProteinSequence(String seqString, CompoundSet<AminoAcidCompound> compoundSet) {
        super(seqString, compoundSet);
    }

    /**
     * A protein sequence where the storage of the sequence is somewhere else. Could be
     * loaded from a large Fasta file or via a Uniprot Proxy reader via Uniprot ID
     * @param proxyLoader
     */
    public ProteinSequence(ProxySequenceReader<AminoAcidCompound> proxyLoader) {
        this(proxyLoader, AminoAcidCompoundSet.getAminoAcidCompoundSet());
    }

    /**
     * A protein sequence where the storage of the sequence is somewhere else with user defined
     * set of amino acids. Could be loaded from a large Fasta file or via a Uniprot Proxy reader
     * via Uniprot ID
     * @param proxyLoader
     */
    public ProteinSequence(ProxySequenceReader<AminoAcidCompound> proxyLoader, CompoundSet<AminoAcidCompound> compoundSet) {
        super(proxyLoader, compoundSet);
    }

    /**
     * A Protein sequence can be stand alone or loaded from a transcript sequence. The design goal is to allow the creation
     * of a Protein sequence from a Uniprot ID or some other Protein ID that based on cross reference you should be able to
     * get the GeneSequence that codes for the protein if the CDS/Gene region is known. From the GeneSequence you should then
     * be able to get the ChromosomeSequence which then allows you explore flaning regions of the gene sequences. The
     * framework is in place to do this but currently hasn't been implement in the reverse direction starting from the
     * Protein sequence.
     *
     * @param parentDNASequence
     * @param begin
     * @param end
     */
    public void setParentDNASequence(AbstractSequence parentDNASequence, Integer begin, Integer end) {
        this.setParentSequence(parentDNASequence);
        setBioBegin(begin);
        setBioEnd(end);
    }

    public static void main(String[] args) {
        ProteinSequence proteinSequence = new ProteinSequence("ARNDCEQGHILKMFPSTWYVBZJX");
        System.out.println(proteinSequence.toString());

        StringProxySequenceReader<AminoAcidCompound> sequenceStringProxyLoader = new StringProxySequenceReader<AminoAcidCompound>("XRNDCEQGHILKMFPSTWYVBZJA", AminoAcidCompoundSet.getAminoAcidCompoundSet());
        ProteinSequence proteinSequenceFromProxy = new ProteinSequence(sequenceStringProxyLoader);
        System.out.println(proteinSequenceFromProxy.toString());

    }
}
