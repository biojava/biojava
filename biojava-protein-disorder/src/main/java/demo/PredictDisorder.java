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
package demo;

import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.loader.UniprotProxySequenceReader;
import org.biojava.nbio.ronn.Jronn;

import java.util.Arrays;

public class PredictDisorder {

	public static void main(String[] args) throws Exception{

		String uniprotID = "O30642";

		ProteinSequence seq = getUniprot(uniprotID);
		System.out.println("Protein Sequence: "+ seq.toString());
		AminoAcidCompoundSet compoundSet = AminoAcidCompoundSet.getAminoAcidCompoundSet();

		if (!compoundSet.isValidSequence(seq) ) {
			System.err.println("Invalid sequence, exiting");
			System.exit(1);
		}

		float[] values = Jronn.getDisorderScores(seq);

		System.out.println("Disorder Scores: "+ Arrays.toString(values));


	}

	/** Fetch a protein sequence from the UniProt web site
	 *
	 * @param uniProtID
	 * @return a Protein Sequence
	 * @throws Exception
	 */
	private static ProteinSequence getUniprot(String uniProtID) throws Exception {

		AminoAcidCompoundSet set = AminoAcidCompoundSet.getAminoAcidCompoundSet();
		UniprotProxySequenceReader<AminoAcidCompound> uniprotSequence = new UniprotProxySequenceReader<AminoAcidCompound>(uniProtID,set);

		ProteinSequence seq = new ProteinSequence(uniprotSequence);

		return seq;
	}
}
