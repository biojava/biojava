package demo;

import java.util.Arrays;

import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.loader.UniprotProxySequenceReader;
import org.biojava3.ronn.Jronn;

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
