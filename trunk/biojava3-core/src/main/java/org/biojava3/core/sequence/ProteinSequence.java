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
 *
 * @author Scooter Willis
 */
public class ProteinSequence extends AbstractSequence<AminoAcidCompound> {



    public ProteinSequence(String seqString) {
        this(seqString, AminoAcidCompoundSet.getAminoAcidCompoundSet());
    }

    public ProteinSequence(String seqString, CompoundSet<AminoAcidCompound> compoundSet) {
        super(seqString, compoundSet);
    }

    public ProteinSequence(ProxySequenceReader<AminoAcidCompound> proxyLoader) {
        this(proxyLoader, AminoAcidCompoundSet.getAminoAcidCompoundSet());
    }

    public ProteinSequence(ProxySequenceReader<AminoAcidCompound> proxyLoader, CompoundSet<AminoAcidCompound> compoundSet) {
        super(proxyLoader, compoundSet);
    }

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
