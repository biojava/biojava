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
 */
package org.biojava.nbio.core.search.io;

import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

/**
 * This class models a search high scoring pair (HSP).
 * You will retrieve a list of this using iterator of a Hit
 *
 * Designed by Paolo Pavan.
 * You may want to find my contacts on Github and LinkedIn for code info
 * or discuss major changes.
 * https://github.com/paolopavan
 *
 * @author Paolo Pavan
 */

public abstract class Hsp <S extends Sequence<C>, C extends Compound> {
	//private static final Logger logger = LoggerFactory.getLogger(Hsp.class);
	private Integer hspNum;
	private Double hspBitScore;
	private Integer hspScore;
	private Double hspEvalue;
	private Integer hspQueryFrom;
	private Integer hspQueryTo;
	private Integer hspHitFrom;
	private Integer hspHitTo;
	private Integer hspQueryFrame;
	private Integer hspHitFrame;
	private Integer hspIdentity;
	private Integer hspPositive;
	private Integer hspGaps;
	private Integer hspAlignLen;
	private Double percentageIdentity = null;
	private Integer mismatchCount = null;
	private SequencePair<S, C> alignment;

	@Override
	public int hashCode() {
		return alignment.hashCode();
	}
	/**
	 * Experimental.
	 * Wants to implement conceptual comparisons of search results.
	 * Fields unrelated to search are deliberately not considered.
	 *
	 * In HSP case, alignment representation strings are considered.
	 * @return true if HSP alignments are the same,
	 * false otherwise or if alignment strings are undetermined
	 */
	@Override
	public boolean equals(Object obj) {
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		final Hsp<?, ?> other = (Hsp<?, ?>) obj;
		return alignment.equals(other.alignment);
	}

	public SequencePair<S,C> getAlignment(){
		return alignment;
	}

//	private static Class<? extends Sequence<?>> guessSequenceType(String gappedSequenceString) {
//		String sequenceString = gappedSequenceString.replace("-", "");
//
//		if (sequenceString.matches("^[ACTG]+$"))
//			return DNASequence.class;
//		else if (sequenceString.matches("^[ACUG]+$"))
//			return RNASequence.class;
//		else
//			return ProteinSequence.class;
//	}
//	private static Sequence<?> getSequence(String gappedSequenceString){
//		if (gappedSequenceString == null) return null;
//
//		Sequence<?> returnSeq = null;
//		String sequenceString = gappedSequenceString.replace("-", "");
//
//		try {
//			if (sequenceString.matches("^[ACTG]+$"))
//				returnSeq = new DNASequence(sequenceString, DNACompoundSet.getDNACompoundSet());
//			else if (sequenceString.matches("^[ACUG]+$"))
//				returnSeq = new RNASequence(sequenceString, DNACompoundSet.getDNACompoundSet());
//			else
//				returnSeq = new ProteinSequence(sequenceString, AminoAcidCompoundSet.getAminoAcidCompoundSet());
//		} catch (CompoundNotFoundException ex) {
//			logger.error("Unexpected error, could not find compound when creating Sequence object from Hsp", ex);
//		}
//		return returnSeq;
//	}

	public int getHspNum() {
		return hspNum;
	}

	public double getHspBitScore() {
		return hspBitScore;
	}

	public int getHspScore() {
		return hspScore;
	}

	public double getHspEvalue() {
		return hspEvalue;
	}

	public int getHspQueryFrom() {
		return hspQueryFrom;
	}

	public int getHspQueryTo() {
		return hspQueryTo;
	}

	public int getHspHitFrom() {
		return hspHitFrom;
	}

	public int getHspHitTo() {
		return hspHitTo;
	}

	public int getHspQueryFrame() {
		return hspQueryFrame;
	}

	public int getHspHitFrame() {
		return hspHitFrame;
	}

	public int getHspIdentity() {
		return hspIdentity;
	}

	public int getHspPositive() {
		return hspPositive;
	}

	public int getHspGaps() {
		return hspGaps;
	}

	public int getHspAlignLen() {
		return hspAlignLen;
	}

	public Double getPercentageIdentity() {
		if (percentageIdentity != null) return percentageIdentity;
		if (hspIdentity!= null && hspAlignLen != null) return (double)hspIdentity/hspAlignLen;
		return null;
	}

	public Integer getMismatchCount() {
		if (mismatchCount != null) return mismatchCount;
		if (hspIdentity!= null && hspAlignLen != null) return hspIdentity-hspAlignLen;
		return null;
	}

	public Hsp(int hspNum, double hspBitScore, int hspScore, double hspEvalue,
			int hspQueryFrom, int hspQueryTo, int hspHitFrom, int hspHitTo,
			int hspQueryFrame, int hspHitFrame, int hspIdentity, int hspPositive,
			int hspGaps, int hspAlignLen, 
			SequencePair<S,C> alignment, 
			Double percentageIdentity, Integer mismatchCount) {
		this.hspNum = hspNum;
		this.hspBitScore = hspBitScore;
		this.hspScore = hspScore;
		this.hspEvalue = hspEvalue;
		this.hspQueryFrom = hspQueryFrom;
		this.hspQueryTo = hspQueryTo;
		this.hspHitFrom = hspHitFrom;
		this.hspHitTo = hspHitTo;
		this.hspQueryFrame = hspQueryFrame;
		this.hspHitFrame = hspHitFrame;
		this.hspIdentity = hspIdentity;
		this.hspPositive = hspPositive;
		this.hspGaps = hspGaps;
		this.hspIdentity = hspAlignLen;
		this.alignment = alignment;
		this.percentageIdentity = percentageIdentity;
		this.mismatchCount = mismatchCount;

		if(this.alignment == null) {
			throw new NullPointerException("Null alignment");
		}
		// sanity check
		if (percentageIdentity != null && (percentageIdentity < 0 || percentageIdentity >1))
			throw new IllegalArgumentException("Percentage identity must be between 0 and 1");

	}

}
