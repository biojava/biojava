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
package org.biojava.nbio.core.sequence.storage;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.Strand;
import org.biojava.nbio.core.sequence.template.*;
import org.biojava.nbio.core.util.Equals;
import org.biojava.nbio.core.util.Hashcoder;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Stores a Sequence as a collection of compounds in an ArrayList
 *
 * @param <C>
 */
public class ArrayListSequenceReader<C extends Compound> implements SequenceReader<C> {

	private CompoundSet<C> compoundSet;
	private ArrayList<C> parsedCompounds = new ArrayList<C>();

	private volatile Integer hashcode = null;

	/**
	 *
	 */
	public ArrayListSequenceReader() {
		//Do nothing
	}

	/**
	 *
	 * @param compounds
	 * @param compoundSet
	 */
	public ArrayListSequenceReader(List<C> compounds, CompoundSet<C> compoundSet) {
		setCompoundSet(compoundSet);
		setContents(compounds);
	}

	/**
	 *
	 * @param sequence
	 * @param compoundSet
	 * @throws CompoundNotFoundException
	 */
	public ArrayListSequenceReader(String sequence, CompoundSet<C> compoundSet) throws CompoundNotFoundException {
		setCompoundSet(compoundSet);
		setContents(sequence);
	}

	/**
	 *
	 * @return
	 */
	@Override
	public String getSequenceAsString() {
		return getSequenceAsString(1, getLength(), Strand.POSITIVE);
	}


	/**
	 *
	 * @param begin
	 * @param end
	 * @param strand
	 * @return
	 */
	public String getSequenceAsString(Integer begin, Integer end, Strand strand) {
		// TODO Optimise/cache.
		SequenceAsStringHelper<C> sequenceAsStringHelper = new SequenceAsStringHelper<C>();
		return sequenceAsStringHelper.getSequenceAsString(this.parsedCompounds, compoundSet, begin, end, strand);
	}

	/**
	 *
	 * @return
	 */
	@Override
	public List<C> getAsList() {
		return this.parsedCompounds;
	}

	/**
	 *
	 * @param position
	 * @return
	 */
	@Override
	public C getCompoundAt(int position) {
		return this.parsedCompounds.get(position - 1);
	}

	/**
	 *
	 * @param compound
	 * @return
	 */
	@Override
	public int getIndexOf(C compound) {
		return this.parsedCompounds.indexOf(compound) + 1;
	}

	/**
	 *
	 * @param compound
	 * @return
	 */
	@Override
	public int getLastIndexOf(C compound) {
		return this.parsedCompounds.lastIndexOf(compound) + 1;
	}

	/**
	 *
	 * @return
	 */
	@Override
	public int getLength() {
		return this.parsedCompounds.size();
	}

	/**
	 *
	 * @return
	 */
	@Override
	public Iterator<C> iterator() {
		return this.parsedCompounds.iterator();
	}

	/**
	 *
	 * @param compoundSet
	 */
	@Override
	public void setCompoundSet(CompoundSet<C> compoundSet) {
		this.compoundSet = compoundSet;
	}

	/**
	 *
	 * @return
	 */
	@Override
	public CompoundSet<C> getCompoundSet() {
		return compoundSet;
	}

	/**
	 *
	 * @param sequence
	 */
	@Override
	public void setContents(String sequence) throws CompoundNotFoundException {
		// Horrendously inefficient - pretty much the way the old BJ did things.
		// TODO Should be optimised.
		this.parsedCompounds.clear();
		hashcode = null;
		int maxCompoundLength = compoundSet.getMaxSingleCompoundStringLength();
		boolean maxCompundLengthEqual1 = true;
		if (maxCompoundLength > 1) {
			maxCompundLengthEqual1 = false;
		}
		int length = sequence.length();
		parsedCompounds.ensureCapacity(length); //get the array size correct
		for (int i = 0; i < length;) {
			String compoundStr = null;
			C compound = null;
			if (maxCompundLengthEqual1) { // trying to save some steps where typically the answer is 1 so avoid complicated for loop
				compoundStr = sequence.substring(i, i + 1);
				compound = compoundSet.getCompoundForString(compoundStr);
			} else {
				for (int compoundStrLength = 1; compound == null && compoundStrLength <= maxCompoundLength; compoundStrLength++) {
					compoundStr = sequence.substring(i, i + compoundStrLength);
					compound = compoundSet.getCompoundForString(compoundStr);
				}
			}
			if (compound == null) {
				throw new CompoundNotFoundException("Cannot find compound for: " + compoundStr);
			} else {
				i += compoundStr.length();
			}
			this.parsedCompounds.add(compound);
		}
		parsedCompounds.trimToSize(); // just in case it increases capacity free up extra memory
	}

	/**
	 *
	 * @param list
	 */
	public void setContents(List<C> list) {
		parsedCompounds.clear();
		for (C c : list) {
			parsedCompounds.add(c);
		}
	}

	/**
	 *
	 * @param bioBegin
	 * @param bioEnd
	 * @return
	 */
	@Override
	public SequenceView<C> getSubSequence(final Integer bioBegin, final Integer bioEnd) {
		return new SequenceProxyView<C>(ArrayListSequenceReader.this, bioBegin, bioEnd);
	}

	/**
	 *
	 * @return
	 */
	@Override
	public AccessionID getAccession() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	/**
	 *
	 * @param compounds
	 * @return
	 */
	@Override
	public int countCompounds(C... compounds) {
		return SequenceMixin.countCompounds(this, compounds);
	}

	/**
	 *
	 * @return
	 */
	@Override
	public SequenceView<C> getInverse() {
		return SequenceMixin.inverse(this);
	}

	@Override
	public int hashCode() {
		if(hashcode == null) {
			int s = Hashcoder.SEED;
			s = Hashcoder.hash(s, parsedCompounds);
			s = Hashcoder.hash(s, compoundSet);
			hashcode = s;
		}
		return hashcode;
	}

	@Override
	@SuppressWarnings("unchecked")
	public boolean equals(Object o) {
		if(Equals.classEqual(this, o)) {
			ArrayListSequenceReader<C> that = (ArrayListSequenceReader<C>)o;
			return  Equals.equal(parsedCompounds, that.parsedCompounds) &&
					Equals.equal(compoundSet, that.compoundSet);
		}
		return false;
	}
}
