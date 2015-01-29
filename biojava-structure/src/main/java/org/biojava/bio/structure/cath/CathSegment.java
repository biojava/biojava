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
 * Date:   2012-6-23
 */

package org.biojava.bio.structure.cath;

import java.io.Serializable;

/**
 *
 * @author Daniel Asarnow
 */
public class CathSegment implements Serializable{

    public static final long serialVersionUID = 1L;

    /**
     * The number of this segment within the domain.
     */
    Integer segmentId;

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

    /**
     * FASTA header.
     */
    String sequenceHeader;

    /**
     * FASTA sequence.
     */
    String sequence;

    public Integer getSegmentId() {
        return segmentId;
    }

    public void setSegmentId(Integer segmentId) {
        this.segmentId = segmentId;
    }

    public String getStart() {
        return start;
    }

    public void setStart(String start) {
        this.start = start;
    }

    public String  getStop() {
        return stop;
    }

    public void setStop(String  stop) {
        this.stop = stop;
    }

    public Integer getLength() {
        return length;
    }

    public void setLength(Integer length) {
        this.length = length;
    }

    public String getSequenceHeader() {
        return sequenceHeader;
    }

    public void setSequenceHeader(String sequenceHeader) {
        this.sequenceHeader = sequenceHeader;
    }

    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

	@Override
	public String toString() {
		return "CathSegment [segmentId=" + segmentId + ", start=" + start
				+ ", stop=" + stop + ", length=" + length + ", sequenceHeader="
				+ sequenceHeader + ", sequence=" + sequence + "]";
	}

   
}
