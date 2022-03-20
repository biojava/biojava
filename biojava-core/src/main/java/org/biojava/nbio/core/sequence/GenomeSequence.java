package org.biojava.nbio.core.sequence;

/**
 * @author Elizabeth James
 */
public class GenomeSequence extends DNASequence{
    @Override
    public int getLength() {
        return Math.abs(this.getBioEnd() - this.getBioBegin()) + 1;
    }
}
