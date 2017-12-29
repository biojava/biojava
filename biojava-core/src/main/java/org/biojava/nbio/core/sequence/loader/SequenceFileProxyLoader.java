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
import org.biojava.nbio.core.sequence.io.template.SequenceParserInterface;
import org.biojava.nbio.core.sequence.storage.SequenceAsStringHelper;
import org.biojava.nbio.core.sequence.template.*;
import org.biojava.nbio.core.util.Equals;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * This class represents the storage container of a sequence stored in a fasta file where
 * the initial parsing of the file we store the offset and length of the sequence. When a call
 * is made to any method that needs sequence data then the file will be opened and the sequence
 * loaded. This class could be improved by using the hints or a some algorithm that indicates
 * the sequence data once loaded should stay loaded. Could keep track of the last time sequence
 * data was loaded and then after X amount of time clear the contents to free up memory.
 *
 *
 * @author Scooter Willis <willishf at gmail dot com>
 * @param <C>
 */
public class SequenceFileProxyLoader<C extends Compound> implements ProxySequenceReader<C> {

	SequenceParserInterface sequenceParser;
	private CompoundSet<C> compoundSet;
	private List<C> parsedCompounds = new ArrayList<C>();
	File file;
	long sequenceStartIndex = -1;
	int sequenceLength = -1;
	//private boolean initialized = false;

	/**
	 *
	 * @param file The file where the sequence will be found
	 * @param sequenceParser The parser to use to load the sequence
	 * @param sequenceStartIndex The file offset to the start of the sequence
	 * @param sequenceLength The length of the sequence
	 * @param compoundSet
	 * @throws IOException if problems occur while reading the file
	 * @throws CompoundNotFoundException if a compound in the sequence can't be found in the given compoundSet
	 */
	public SequenceFileProxyLoader(File file, SequenceParserInterface sequenceParser, long sequenceStartIndex, int sequenceLength, CompoundSet<C> compoundSet)
			throws IOException, CompoundNotFoundException {
		this.sequenceParser = sequenceParser;
		this.file = file;
		this.sequenceStartIndex = sequenceStartIndex;
		this.sequenceLength = sequenceLength;
		setCompoundSet(compoundSet);

		init();
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
	 *  Load the sequence
	 * @return
	 */
	private boolean init() throws IOException, CompoundNotFoundException {

		BufferedReader br = new BufferedReader(new FileReader(file));
		br.skip(sequenceStartIndex);
		String sequence = sequenceParser.getSequence(br, sequenceLength);
		setContents(sequence);
		br.close(); // close file to prevent too many being open

		return true;
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

	/**
	 *
	 * @return
	 */
	@Override
	public int getLength() {
		return sequenceLength;
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
	public String toString() {

		return getSequenceAsString();
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
	 * @param bioBegin
	 * @param bioEnd
	 * @param strand
	 * @return
	 */
	public String getSequenceAsString(Integer bioBegin, Integer bioEnd, Strand strand) {

		SequenceAsStringHelper<C> sequenceAsStringHelper = new SequenceAsStringHelper<C>();
		return sequenceAsStringHelper.getSequenceAsString(this.parsedCompounds, compoundSet, bioBegin, bioEnd, strand);
	}

	/**
	 *
	 * @return
	 */
	@Override
	public List<C> getAsList() {

		return this.parsedCompounds;

	}

	@Override
	public boolean equals(Object o) {

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

	/**
	 *
	 * @param bioBegin
	 * @param bioEnd
	 * @return
	 */
	@Override
	public SequenceView<C> getSubSequence(final Integer bioBegin, final Integer bioEnd) {

		return new SequenceProxyView<C>(SequenceFileProxyLoader.this, bioBegin, bioEnd);
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
	 * @return
	 */
	@Override
	public CompoundSet<C> getCompoundSet() {
		return compoundSet;
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
}
