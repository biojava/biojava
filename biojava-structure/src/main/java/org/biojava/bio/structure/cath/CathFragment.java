/*
 * BioJava development code
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
 * Author: Daniel Asarnow
 * Date:   2012-7-23
 */

package org.biojava.bio.structure.cath;

import java.io.Serializable;

/**
 * @author Daniel Asarnow
 */
public class CathFragment implements Serializable{

    public static final long serialVersionUID = 1L;

    /**
     * The number of this segment within the domain.
     */
    Integer fragmentId;

    /**
     * The first residue in the segment.
     * Refers to the complete residue specification (sequence number AND insertion code).
     */
    String start;

     /**
     * The last residue in the segment.
     * Refers to the complete residue specification (sequence number AND insertion code).
     */
    String stop;

    /**
     * Number of residues in the segment. This value is parsed, not calculated.
     */
    Integer length;

    public Integer getFragmentId() {
        return fragmentId;
    }

    public void setFragmentId(Integer fragmentId) {
        this.fragmentId = fragmentId;
    }

    public String getStart() {
        return start;
    }

    public void setStart(String start) {
        this.start = start;
    }

    public String getStop() {
        return stop;
    }

    public void setStop(String stop) {
        this.stop = stop;
    }

    public Integer getLength() {
        return length;
    }

    public void setLength(Integer length) {
        this.length = length;
    }
}
