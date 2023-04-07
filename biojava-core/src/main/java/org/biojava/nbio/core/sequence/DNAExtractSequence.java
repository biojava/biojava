package org.biojava.nbio.core.sequence;

public class DNAExtractSequence extends DNASequence{
    @Override
    public int getLength() {
        return Math.abs(this.getBioEnd() - this.getBioBegin()) + 1;
    }
}
