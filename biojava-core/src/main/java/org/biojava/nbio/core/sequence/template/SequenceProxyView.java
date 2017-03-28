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
 *
 * @author Richard Holland
 *
 *
 */
package org.biojava.nbio.core.sequence.template;

import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.util.Equals;

import java.util.Iterator;
import java.util.List;

public class SequenceProxyView<C extends Compound> implements SequenceView<C> {

	private Integer bioStart;
	private Integer bioEnd;
	private Sequence<C> sequence;

	public SequenceProxyView() {
	}

	public SequenceProxyView(Sequence<C> sequence) {
		this(sequence, 1, sequence.getLength());
	}

	/**
	 * Main constructor for working with SequenceProxyViews
	 *
	 * @param sequence Sequence to proxy
	 * @param bioStart Start; cannot be less than 1
	 * @param bioEnd End; cannot be greater than the sequence length
	 */
	public SequenceProxyView(Sequence<C> sequence, Integer bioStart, Integer bioEnd) {
		this.sequence = sequence;
		setBioStart(bioStart);
		setBioEnd(bioEnd);
	}

	@Override
	public Sequence<C> getViewedSequence() {
		return sequence;
	}

	@Override
	public String getSequenceAsString() {
		return SequenceMixin.toString(this);
	}

	@Override
	public List<C> getAsList() {
		return SequenceMixin.toList(this);
	}

	@Override
	public boolean equals(Object o) {

		if(! Equals.classEqual(this, o)) {
			return false;
		}

		Sequence<C> other = (Sequence<C>)o;

		return sequence.equals(other);
	}

	@Override
	public int hashCode(){
		String s = getSequenceAsString();
		return s.hashCode();
	}

	@Override
	public C getCompoundAt(int position) {
		return getViewedSequence().getCompoundAt((getBioStart() + position) - 1);
	}

	@Override
	public int getIndexOf(C compound) {
		return SequenceMixin.indexOf(this, compound);
	}

	@Override
	public int getLastIndexOf(C compound) {
		return SequenceMixin.lastIndexOf(this, compound);
	}

	@Override
	public int getLength() {
		return (getBioEnd() - getBioStart()) + 1;
	}

	@Override
	public CompoundSet<C> getCompoundSet() {
		return getViewedSequence().getCompoundSet();
	}

	@Override
	public SequenceView<C> getSubSequence(final Integer bioStart, final Integer bioEnd) {
		return new SequenceProxyView<C>(this, bioStart, bioEnd);
	}

	@Override
	public Iterator<C> iterator() {
		return new SequenceMixin.SequenceIterator<C>(this);
	}

	@Override
	public AccessionID getAccession() {
		return getViewedSequence().getAccession();
	}

	/**
	 * @return the bioStart
	 */
	@Override
	public Integer getBioStart() {
		return bioStart;
	}

	/**
	 * @param bioStart the bioStart to set
	 */
	public void setBioStart(Integer bioStart) {
		if (bioStart < 1) {
			throw new IllegalArgumentException("The given start "
					+ bioStart + " is less than 1; cannot index less than 1");
		}
		this.bioStart = bioStart;
	}

	/**
	 * @return the bioEnd
	 */
	@Override
	public Integer getBioEnd() {
		return bioEnd;
	}

	/**
	 * @param bioEnd the bioEnd to set
	 */
	public void setBioEnd(Integer bioEnd) {
		if(sequence == null) {
			throw new NullPointerException("No sequence given before setting the end coordinate; cannot be done");
		}
		// had a bug in the code that allowed this to work. The length of a any exon or cds sequence was always the length of the
		//parent sequence. Sequence class doesn't have bioStart and bioEnd exposed to do a proper comparison of getting
		// a subsequence. Januar-20=2011 Scooter
	 //   if (bioEnd > sequence.getLength()) {
	 //       throw new IllegalArgumentException("The given end "
	 //               + bioEnd + " is greater than sequence length"
	 //               + sequence.getLength());
	 //   }
		this.bioEnd = bioEnd;
	}

	@Override
	public int countCompounds(C... compounds) {
		return SequenceMixin.countCompounds(this, compounds);
	}

	@Override
	public SequenceView<C> getInverse() {
		return SequenceMixin.inverse(this);
	}
}
