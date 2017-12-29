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
 * @auther Scooter Willis
 *
 */
package org.biojava.nbio.core.sequence.loader;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.Strand;
import org.biojava.nbio.core.sequence.storage.SequenceAsStringHelper;
import org.biojava.nbio.core.sequence.template.*;
import org.biojava.nbio.core.util.Equals;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * An example of a ProxySequenceReader that is created from a String. Used for testing
 * @author Scooter Willis <willishf at gmail dot com>
 * @param <C>
 */

public class StringProxySequenceReader<C extends Compound> implements ProxySequenceReader<C> {

	private String sequence;
	private CompoundSet<C> compoundSet;
	private List<C> parsedCompounds = new ArrayList<C>();

	public StringProxySequenceReader() {}

	public StringProxySequenceReader(String sequence, CompoundSet<C> compoundSet) throws CompoundNotFoundException {
		this.sequence = sequence;
		setCompoundSet(compoundSet);
		setContents(sequence);
	}

	@Override
	public void setCompoundSet(CompoundSet<C> compoundSet) {
		this.compoundSet = compoundSet;
	}

	@Override
	public void setContents(String sequence) throws CompoundNotFoundException {
		// Horrendously inefficient - pretty much the way the old BJ did things.
		// TODO Should be optimised.
		this.sequence = sequence;
		this.parsedCompounds.clear();
		for (int i = 0; i < sequence.length();) {
			String compoundStr = null;
			C compound = null;
			for (int compoundStrLength = 1; compound == null && compoundStrLength <= compoundSet.getMaxSingleCompoundStringLength(); compoundStrLength++) {
				compoundStr = sequence.substring(i, i + compoundStrLength);
				compound = compoundSet.getCompoundForString(compoundStr);
			}
			if (compound == null) {
				throw new CompoundNotFoundException("Compound "+compoundStr+" not found");
			} else {
				i += compoundStr.length();
			}
			this.parsedCompounds.add(compound);
		}
	}

	public void setContents(String sequence, ArrayList features) throws CompoundNotFoundException{
		setContents(sequence);
	}

	@Override
	public int getLength() {
		return this.parsedCompounds.size();
	}

	@Override
	public C getCompoundAt(int position) {
		return this.parsedCompounds.get(position - 1);
	}

	@Override
	public int getIndexOf(C compound) {
		return this.parsedCompounds.indexOf(compound) + 1;
	}

	@Override
	public int getLastIndexOf(C compound) {
		return this.parsedCompounds.lastIndexOf(compound) + 1;
	}


	@Override
	public String toString() {
		return getSequenceAsString();
	}

	@Override
	public String getSequenceAsString() {
		return sequence;
	}

	@Override
	public List<C> getAsList() {
		return this.parsedCompounds;
	}



	public String getSequenceAsString(Integer bioBegin, Integer bioEnd,Strand strand) {
		SequenceAsStringHelper<C> sequenceAsStringHelper = new SequenceAsStringHelper<C>();
		return sequenceAsStringHelper.getSequenceAsString(this.parsedCompounds, compoundSet, bioBegin, bioEnd, strand);
	}

	@Override
	public SequenceView<C> getSubSequence(final Integer bioBegin, final Integer bioEnd) {
		return new SequenceProxyView<C>(StringProxySequenceReader.this,bioBegin,bioEnd);
	}

	@Override
	public Iterator<C> iterator() {
		return this.parsedCompounds.iterator();
	}

	@Override
	public CompoundSet<C> getCompoundSet() {
	  return compoundSet;
	}


	@Override
	public AccessionID getAccession() {
		throw new UnsupportedOperationException("Not supported yet.");
	}


	@Override
	public int countCompounds(C... compounds) {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public SequenceView<C> getInverse() {
		return SequenceMixin.inverse(this);
	}

	@Override
	public boolean equals(Object o){

		if(! Equals.classEqual(this, o)) {
			return false;
		}

		Sequence<C> other = (Sequence<C>)o;
		if ( other.getCompoundSet() != getCompoundSet())
			return false;

		List<C> rawCompounds = getAsList();
		List<C> otherCompounds = other.getAsList();

		if ( rawCompounds.size() != otherCompounds.size())
			return false;

		for (int i = 0 ; i < rawCompounds.size() ; i++){
			Compound myCompound = rawCompounds.get(i);
			Compound otherCompound = otherCompounds.get(i);
			if ( ! myCompound.equalsIgnoreCase(otherCompound))
				return false;
		}
		return true;
	}

	@Override
	public int hashCode(){
		String s = getSequenceAsString();
		return s.hashCode();
	}
}
