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
 * Created on November 21, 2010
 * Author: Mark Chapman
 */

package org.biojava.nbio.core.sequence;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.LightweightProfile;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.EnumSource;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import static org.junit.jupiter.api.Assertions.*;


class MultipleSequenceAlignmentTest {

	private MultipleSequenceAlignment<ProteinSequence, AminoAcidCompound> msaProteins;
	private MultipleSequenceAlignment<DNASequence,NucleotideCompound> msaDNA;

	private static final String aaSeq = "ARNDCEQGHILKMFPSTWYVBZJX";
	@BeforeEach
	 void setup() throws CompoundNotFoundException {
		msaProteins = new MultipleSequenceAlignment<>();
		for (int i = 0; i < 8; i++) {
			ProteinSequence ps = new ProteinSequence(aaSeq);
			ps.setAccession(new AccessionID(i+""));
			msaProteins.addAlignedSequence(ps);
		}
		msaDNA = new MultipleSequenceAlignment<>();
		for (int i = 0; i < 7; i++) {
			msaDNA.addAlignedSequence(new DNASequence("ATCGATCGATCGATCG"));
		}
	}

	@Test
	void allSequencesMustBeSameLength() throws CompoundNotFoundException {
		ProteinSequence differentLength = new ProteinSequence("ARNDC");
		assertThrows(IllegalArgumentException.class, ()->msaProteins.addAlignedSequence(differentLength));
	}

	@Test
	void addRemoveAlignments() throws CompoundNotFoundException {
		assertEquals(8, msaProteins.getSize());
		assertEquals(8, msaProteins.getAlignedSequences().size());
		assertEquals(aaSeq.length(), msaProteins.getLength());
		msaProteins.removeAlignedSequence(new ProteinSequence(aaSeq));
		assertEquals(7, msaProteins.getSize());
		assertEquals(7, msaProteins.getAlignedSequences().size());
	}

	@ParameterizedTest
	@EnumSource(LightweightProfile.StringFormat.class)
	void formattedAlignmentToString(LightweightProfile.StringFormat format){
		String formatted = msaProteins.toString(format);
		assertTrue(formatted.length() > 0);
	}

	@Test
	void alignmentToBasicString(){
		String alnStr = msaProteins.toString();
		String [] lines  = alnStr.split(System.lineSeparator());
		assertEquals(8, lines.length);

		//lines all same length
		Set<Integer> collect = Arrays.stream(lines).map(String::length).collect(Collectors.toSet());
		assertEquals(1, collect.size());
	}
	@Test
	void alignmentToWidth() {
		String alnStr = msaProteins.toString(10);
		assertEquals(29, alnStr.split(System.lineSeparator()).length);
	}

	@Test
	 void testGetCompoundsAt() {
		AminoAcidCompound aminoAcid = AminoAcidCompoundSet.getAminoAcidCompoundSet().getCompoundForString("N");
		List<AminoAcidCompound> colProteins = new ArrayList<>();
		for (int i = 0; i < 8; i++) {
			colProteins.add(aminoAcid);
		}
		assertEquals(msaProteins.getCompoundsAt(3), colProteins);
		NucleotideCompound nucleotide = DNACompoundSet.getDNACompoundSet().getCompoundForString("C");
		List<NucleotideCompound> colDNA = new ArrayList<>();
		for (int i = 0; i < 7; i++) {
			colDNA.add(nucleotide);
		}
		assertEquals(msaDNA.getCompoundsAt(3), colDNA);
	}


}
