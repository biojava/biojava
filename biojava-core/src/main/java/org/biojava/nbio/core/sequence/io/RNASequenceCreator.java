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

package org.biojava.nbio.core.sequence.io;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.RNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava.nbio.core.sequence.loader.ArrayListProxySequenceReader;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.ProxySequenceReader;

import java.util.List;

/**
 * Used to create a RNA sequence
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class RNASequenceCreator implements
		SequenceCreatorInterface<NucleotideCompound> {

	private final CompoundSet<NucleotideCompound> compoundSet;
/**
 *
 * @param compoundSet
 */
	public RNASequenceCreator(CompoundSet<NucleotideCompound> compoundSet) {
		this.compoundSet = compoundSet;
	}
/**
 *
 * @param sequence
 * @param index
 * @return
 * @throws CompoundNotFoundException
 */
	@Override
public AbstractSequence<NucleotideCompound> getSequence(String sequence, long index) throws CompoundNotFoundException {
		return new RNASequence(sequence, compoundSet);
	}
/**
 *
 * @param proxyLoader
 * @param index
 * @return
 */
	@Override
public AbstractSequence<NucleotideCompound> getSequence(
			ProxySequenceReader<NucleotideCompound> proxyLoader, long index) {
		return new RNASequence(proxyLoader, compoundSet);
	}
/**
 *
 * @param list
 * @return
 */
	@Override
public AbstractSequence<NucleotideCompound> getSequence(List<NucleotideCompound> list) {
		ArrayListProxySequenceReader<NucleotideCompound> store =
			new ArrayListProxySequenceReader<NucleotideCompound>();
		store.setCompoundSet(compoundSet);
		store.setContents(list);
		return new RNASequence(store);
	}
}
