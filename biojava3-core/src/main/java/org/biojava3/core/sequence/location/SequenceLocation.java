/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence.location;

import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class SequenceLocation<S extends AbstractSequence<C>, C extends Compound> extends SimpleLocation {
private S sequence;
    public SequenceLocation(int start, int end,S sequence){
        super(start,end);
        this.sequence = sequence;

    }

    /**
     * @return the sequence
     */
    public S getSequence() {
        return sequence;
    }


}
