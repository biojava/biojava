/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence.io;

import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.SequenceProxyLoader;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class DNASequenceCreator implements SequenceCreatorInterface{
    CompoundSet compoundSet = null;
    public DNASequenceCreator(CompoundSet compoundSet){
        this.compoundSet = compoundSet;
    }

    public AbstractSequence getSequence(String sequence,long index) {
        return new DNASequence(sequence,compoundSet);
    }

    public AbstractSequence getSequence(SequenceProxyLoader proxyLoader,long index) {
        return new DNASequence(proxyLoader,compoundSet);
    }

}
