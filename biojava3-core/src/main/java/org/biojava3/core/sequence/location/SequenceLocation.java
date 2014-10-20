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

import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.location.template.Location;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;
import java.util.List;
/**
 * A location in a sequence that keeps a reference to its parent sequence
 * @author Scooter Willis <willishf at gmail dot com>
 * @author Paolo Pavan
 */
public class SequenceLocation<S extends AbstractSequence<C>, C extends Compound> extends SimpleLocation {
    private S sequence;

    private boolean partialOn5prime = false;
    private boolean partialOn3prime = false;
    
    public SequenceLocation(int start, int end,S sequence){
        super(start,end);
        this.sequence = sequence;

    }


    public SequenceLocation(int start, int end, S sequence, Strand strand, boolean circular, List<Location> subLocations) {
        super(new SimplePoint(start), new SimplePoint(end), strand, circular, subLocations);

        this.sequence = sequence;
    }

    public SequenceLocation(int start, int end,S sequence, Strand strand){
        super(start,end);
        this.sequence = sequence;
        setStrand(strand);

    }

    /**
     * @return the sequence
     */
    public S getSequence() {
        return sequence;
    }
    
    public void setSequence(S sequence) {
        this.sequence = sequence;
    }
    
    
    public boolean isPartialOn5prime() {
        return partialOn5prime;
    }

    public void setPartialOn5prime(boolean partialOn5prime) {
        this.partialOn5prime = partialOn5prime;
    }

    public boolean isPartialOn3prime() {
        return partialOn3prime;
    }

    public void setPartialOn3prime(boolean partialOn3prime) {
        this.partialOn3prime = partialOn3prime;
    }
    
    public boolean isPartial() {
        return partialOn5prime || partialOn3prime;
    }


}
