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
package org.biojava.nbio.ronn;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.loader.UniprotProxySequenceReader;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

/**
 * Tests the following UniProt IDs:
 *
 * Methyltransferase mtmB1, UniProt O30642
 * Sec (U), Selenoprotein V, UniProt P59797
 * Asx (B), 5-hydroxytryptamine receptor 1F, UniProt Q9N2D9
 * Xaa (X), DNA-3-methyladenine glycosylase (Fragment), UniProt P23571
 * Glx (Z), 14-3-3 protein gamma-1, UniProt Q6UFZ3
 *
 *
 */
public class NonstandardProteinCompoundTest  {

	@Test
	public void test() throws IOException, InterruptedException {


		String[] uniprotIDs = new String[]{"O30642","P59797","Q9N2D9","Q9N2D9","Q6UFZ3","O30642"};

		for ( String id: uniprotIDs){
			try {
				testUniprot(id);
				// throttle load on uniprot server
				Thread.sleep(1000);

			} catch (CompoundNotFoundException e){

				fail(e.getMessage());

			}
		}
	}

	private void testUniprot(String uniprotID) throws CompoundNotFoundException, IOException {

		ProteinSequence seq = getUniprot(uniprotID);

		AminoAcidCompoundSet compoundSet = AminoAcidCompoundSet.getAminoAcidCompoundSet();

/*		for (AminoAcidCompound compound : seq) {
			System.out.println(compound.getShortName() + " " + compound.getLongName() + " " + compound.getDescription() + " | " + compoundSet.getEquivalentCompounds(compound) + " " + compound.getMolecularWeight() + " " + compound.getBase());
		}
		*/
		assertTrue(compoundSet.isValidSequence(seq));



		Jronn.getDisorderScores(seq);


	}

	/**
	 * Fetch a protein sequence from the UniProt web site
	 *
	 * @param uniProtID
	 * @return a Protein Sequence
	 * @throws IOException
	 * @throws CompoundNotFoundException
	 */
	private static ProteinSequence getUniprot(String uniProtID) throws CompoundNotFoundException, IOException {

		AminoAcidCompoundSet set = AminoAcidCompoundSet.getAminoAcidCompoundSet();
		UniprotProxySequenceReader<AminoAcidCompound> uniprotSequence =
				new UniprotProxySequenceReader<AminoAcidCompound>(uniProtID,set);

		ProteinSequence seq = new ProteinSequence(uniprotSequence);

		return seq;
	}

}
