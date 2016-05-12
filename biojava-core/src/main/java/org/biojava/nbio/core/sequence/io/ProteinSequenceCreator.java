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
 */

package org.biojava.nbio.core.sequence.io;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava.nbio.core.sequence.loader.ArrayListProxySequenceReader;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.ProxySequenceReader;

import java.util.List;

/**
 * Used to create a ProteinSequence from a String to allow for details
 * about the location of the sequence etc.
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class ProteinSequenceCreator implements
		SequenceCreatorInterface<AminoAcidCompound> {

	private CompoundSet<AminoAcidCompound> compoundSet;
/**
 *
 * @param compoundSet
 */
	public ProteinSequenceCreator(CompoundSet<AminoAcidCompound> compoundSet) {
		this.compoundSet = compoundSet;
	}
/**
 *
 * @param sequence
 * @param index not used in this implementation
 * @return
 * @throws CompoundNotFoundException
 */
	@Override
public AbstractSequence<AminoAcidCompound> getSequence(String sequence,
			long index) throws CompoundNotFoundException {
		return new ProteinSequence(sequence, compoundSet);
	}
/**
 *
 * @param list
 * @return
 */
	@Override
public AbstractSequence<AminoAcidCompound> getSequence(
			List<AminoAcidCompound> list) {
		ArrayListProxySequenceReader<AminoAcidCompound> store = new ArrayListProxySequenceReader<AminoAcidCompound>();
		store.setCompoundSet(compoundSet);
		store.setContents(list);
		return new ProteinSequence(store);
	}
/**
 *
 * @param proxyLoader
 * @param index not used in this implementation
 * @return
 */
	@Override
public AbstractSequence<AminoAcidCompound> getSequence(
			ProxySequenceReader<AminoAcidCompound> proxyLoader, long index) {
		return new ProteinSequence(proxyLoader, compoundSet);
	}
}
