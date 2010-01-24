/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence.io.template;

import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.SequenceProxyLoader;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public interface SequenceCreatorInterface<C extends Compound> {

    public AbstractSequence<C> getSequence(String sequence, long index);
    public AbstractSequence<C> getSequence(SequenceProxyLoader<C> proxyLoader,long index);
}
