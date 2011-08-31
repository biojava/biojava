/*  @(#)FastaSequence.java 1.0  September 2009
 * 
 *  Copyright (c) 2009 Peter Troshin
 *  
 *        BioJava development code
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
 */

package org.biojava3.data.sequence;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;

/**
 * A FASTA formatted sequence. Please note that this class does not make any
 * assumptions as to what sequence it store e.g. it could be nucleotide, protein
 * or even gapped alignment sequence! The only guarantee it makes is that the
 * sequence does not contain white space characters e.g. spaces, new lines etc
 * 
 * @author pvtroshin
 * @version 1.0
 * @since 3.0.2
 */

@XmlAccessorType(XmlAccessType.FIELD)
public final class FastaSequence implements Comparable<FastaSequence>{

    /**
     * Sequence id
     */
    private String id;

    // TODO what about gapped sequence here! should be indicated
    /**
     * Returns the string representation of sequence
     */
    private String sequence;

    private FastaSequence() {
	// Default constructor for JaxB
    }

    /**
     * Upon construction the any whitespace characters are removed from the
     * sequence
     * 
     * @param id
     * @param sequence
     */
    public FastaSequence(final String id, final String sequence) {
	this.id = id.trim();
	this.sequence = SequenceUtil.cleanSequence(sequence);
    }

    /**
     * Gets the value of id
     * 
     * @return the value of id
     */
    public String getId() {
	return id;
    }

    /**
     * Gets the value of sequence
     * 
     * @return the value of sequence
     */
    public String getSequence() {
	return sequence;
    }

    public static int countMatchesInSequence(final String theString,
	    final String theRegExp) {
	final Pattern p = Pattern.compile(theRegExp);
	final Matcher m = p.matcher(theString);
	int cnt = 0;
	while (m.find()) {
	    cnt++;
	}
	return cnt;
    }

    public String getFormattedFasta() {
	return getFormatedSequence(80);
    }

    /**
     * 
     * @return one line name, next line sequence, no matter what the sequence
     *         length is
     */
    public String getOnelineFasta() {
	String fasta = ">" + getId() + "\n";
	fasta += getSequence() + "\n";
	return fasta;
    }

    /**
     * Format sequence per width letter in one string. Without spaces.
     * 
     * @return multiple line formated sequence, one line width letters length
     * 
     */
    public String getFormatedSequence(final int width) {
	if (sequence == null) {
	    return "";
	}

	assert width >= 0 : "Wrong width parameter ";

	final StringBuilder sb = new StringBuilder(sequence);
	int nchunks = sequence.length() / width;
	// add up inserted new line chars
	nchunks = (nchunks + sequence.length()) / width;
	int nlineCharcounter = 0;
	for (int i = 1; i <= nchunks; i++) {
	    final int insPos = width * i + nlineCharcounter;
	    // to prevent inserting new line in the very end of a sequence then
	    // it would have failed.
	    // Also covers the case when the sequences shorter than width
	    if (sb.length() <= insPos) {
		break;
	    }
	    sb.insert(insPos, "\n");
	    nlineCharcounter++;
	}
	return sb.toString();
    }

    /**
     * 
     * @return sequence length
     */
    public int getLength() {
	return sequence.length();
    }

    /**
     * Same as oneLineFasta
     */
    @Override
    public String toString() {
	return this.getOnelineFasta();
    }

    @Override
    public int hashCode() {
	final int prime = 31;
	int result = 1;
	result = prime * result + ((id == null) ? 0 : id.hashCode());
	result = prime * result
		+ ((sequence == null) ? 0 : sequence.hashCode());
	return result;
    }

    @Override
    public boolean equals(final Object obj) {
	if (this == obj) {
	    return true;
	}
	if (obj == null) {
	    return false;
	}
	if (getClass() != obj.getClass()) {
	    return false;
	}
	final FastaSequence other = (FastaSequence) obj;
	if (id == null) {
	    if (other.id != null) {
		return false;
	    }
	} else if (!id.equals(other.id)) {
	    return false;
	}
	if (sequence == null) {
	    if (other.sequence != null) {
		return false;
	    }
	} else if (!sequence.equals(other.sequence)) {
	    return false;
	}
	return true;
    }

	@Override
	public int compareTo(FastaSequence o) {
		if(o==null || o.id==null) 
			return 1;
		
		return this.getId().compareTo(o.id);
	}

}
