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

/**
 * @author Scooter Willis
 *
 */
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
