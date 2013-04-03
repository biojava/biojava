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
 * Created on 01-21-2010
 */

package org.biojava3.core.sequence.location;

import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;

/**
 * A location in a sequence that keeps a reference to its parent sequence
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
