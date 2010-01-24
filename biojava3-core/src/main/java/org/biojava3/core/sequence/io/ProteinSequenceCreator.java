/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence.io;

import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.SequenceProxyLoader;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class ProteinSequenceCreator implements SequenceCreatorInterface{
    CompoundSet compoundSet = null;
    public ProteinSequenceCreator(CompoundSet compoundSet){
        this.compoundSet = compoundSet;
    }

    public AbstractSequence getSequence(String sequence,long index) {
        return new ProteinSequence(sequence,compoundSet);
    }

    public AbstractSequence getSequence(SequenceProxyLoader proxyLoader,long index) {
        return new ProteinSequence(proxyLoader,compoundSet);
    }

}
