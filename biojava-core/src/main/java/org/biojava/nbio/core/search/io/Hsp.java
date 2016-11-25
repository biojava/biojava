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
		final int prime = 31;
		int result = 1;
		result = prime * result + ((alignment == null) ? 0 : alignment.hashCode());
		result = prime * result + ((hspAlignLen == null) ? 0 : hspAlignLen.hashCode());
		result = prime * result + ((hspBitScore == null) ? 0 : hspBitScore.hashCode());
		result = prime * result + ((hspEvalue == null) ? 0 : hspEvalue.hashCode());
		result = prime * result + ((hspGaps == null) ? 0 : hspGaps.hashCode());
		result = prime * result + ((hspHitFrame == null) ? 0 : hspHitFrame.hashCode());
		result = prime * result + ((hspHitFrom == null) ? 0 : hspHitFrom.hashCode());
		result = prime * result + ((hspHitTo == null) ? 0 : hspHitTo.hashCode());
		result = prime * result + ((hspIdentity == null) ? 0 : hspIdentity.hashCode());
		result = prime * result + ((hspNum == null) ? 0 : hspNum.hashCode());
		result = prime * result + ((hspPositive == null) ? 0 : hspPositive.hashCode());
		result = prime * result + ((hspQueryFrame == null) ? 0 : hspQueryFrame.hashCode());
		result = prime * result + ((hspQueryFrom == null) ? 0 : hspQueryFrom.hashCode());
		result = prime * result + ((hspQueryTo == null) ? 0 : hspQueryTo.hashCode());
		result = prime * result + ((hspScore == null) ? 0 : hspScore.hashCode());
		result = prime * result + ((mismatchCount == null) ? 0 : mismatchCount.hashCode());
		result = prime * result + ((percentageIdentity == null) ? 0 : percentageIdentity.hashCode());
		return result;
	}


	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Hsp<?,?> other = (Hsp<?,?>) obj;
		if (alignment == null) {
			if (other.alignment != null)
				return false;
		} else if (!alignment.equals(other.alignment))
			return false;
		if (hspAlignLen == null) {
			if (other.hspAlignLen != null)
				return false;
		} else if (!hspAlignLen.equals(other.hspAlignLen))
			return false;
		if (hspBitScore == null) {
			if (other.hspBitScore != null)
				return false;
		} else if (!hspBitScore.equals(other.hspBitScore))
			return false;
		if (hspEvalue == null) {
			if (other.hspEvalue != null)
				return false;
		} else if (!hspEvalue.equals(other.hspEvalue))
			return false;
		if (hspGaps == null) {
			if (other.hspGaps != null)
				return false;
		} else if (!hspGaps.equals(other.hspGaps))
			return false;
		if (hspHitFrame == null) {
			if (other.hspHitFrame != null)
				return false;
		} else if (!hspHitFrame.equals(other.hspHitFrame))
			return false;
		if (hspHitFrom == null) {
			if (other.hspHitFrom != null)
				return false;
		} else if (!hspHitFrom.equals(other.hspHitFrom))
			return false;
		if (hspHitTo == null) {
			if (other.hspHitTo != null)
				return false;
		} else if (!hspHitTo.equals(other.hspHitTo))
			return false;
		if (hspIdentity == null) {
			if (other.hspIdentity != null)
				return false;
		} else if (!hspIdentity.equals(other.hspIdentity))
			return false;
		if (hspNum == null) {
			if (other.hspNum != null)
				return false;
		} else if (!hspNum.equals(other.hspNum))
			return false;
		if (hspPositive == null) {
			if (other.hspPositive != null)
				return false;
		} else if (!hspPositive.equals(other.hspPositive))
			return false;
		if (hspQueryFrame == null) {
			if (other.hspQueryFrame != null)
				return false;
		} else if (!hspQueryFrame.equals(other.hspQueryFrame))
			return false;
		if (hspQueryFrom == null) {
			if (other.hspQueryFrom != null)
				return false;
		} else if (!hspQueryFrom.equals(other.hspQueryFrom))
			return false;
		if (hspQueryTo == null) {
			if (other.hspQueryTo != null)
				return false;
		} else if (!hspQueryTo.equals(other.hspQueryTo))
			return false;
		if (hspScore == null) {
			if (other.hspScore != null)
				return false;
		} else if (!hspScore.equals(other.hspScore))
			return false;
		if (mismatchCount == null) {
			if (other.mismatchCount != null)
				return false;
		} else if (!mismatchCount.equals(other.mismatchCount))
			return false;
		if (percentageIdentity == null) {
			if (other.percentageIdentity != null)
				return false;
		} else if (!percentageIdentity.equals(other.percentageIdentity))
			return false;
		return true;
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

		// sanity check
		if (percentageIdentity != null && (percentageIdentity < 0 || percentageIdentity >1))
			throw new IllegalArgumentException("Percentage identity must be between 0 and 1");

	}

}
