package demo;

import java.util.Arrays;

import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.loader.UniprotProxySequenceReader;
import org.biojava3.ronn.Jronn;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class PredictDisorder {
  
	private static final Logger logger = LoggerFactory.getLogger(PredictDisorder.class);

	public static void main(String[] args){

		String uniprotID = "O30642";
		try {
			ProteinSequence seq = getUniprot(uniprotID);
			logger.info("Protein Sequence: {}", seq.toString());
			AminoAcidCompoundSet compoundSet = AminoAcidCompoundSet.getAminoAcidCompoundSet();
			
			compoundSet.verifySequence(seq);
			
			float[] values = Jronn.getDisorderScores(seq);
			
			logger.info("Disorder Scores: {}", Arrays.toString(values));
			
		} catch (Exception e){
			logger.error("Exception: ", e);
		}
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
