package org.biojava3.core.sequence;

import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.compound.RNACompoundSet;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.SequenceProxyLoader;

public class RNASequence extends AbstractSequence<NucleotideCompound> {


    public RNASequence(String seqString) {
        super(seqString, RNACompoundSet.getRNACompoundSet());
    }

    public RNASequence(SequenceProxyLoader<NucleotideCompound> proxyLoader) {
        super(proxyLoader, RNACompoundSet.getRNACompoundSet());
    }


    public RNASequence getReverseComplement() {
        
        
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public RNASequence getReverse() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public RNASequence getComplement() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public DNASequence getDNASequence() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
