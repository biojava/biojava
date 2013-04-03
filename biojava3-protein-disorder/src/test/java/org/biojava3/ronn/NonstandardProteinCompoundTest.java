package org.biojava3.ronn;

import static org.junit.Assert.*;

import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.loader.UniprotProxySequenceReader;
import org.junit.Test;

/**
 * Tests the following UniProt IDs:
 * 
Methyltransferase mtmB1, UniProt O30642
Sec (U), Selenoprotein V, UniProt P59797
Asx (B), 5-hydroxytryptamine receptor 1F, UniProt Q9N2D9
Xaa (X), DNA-3-methyladenine glycosylase (Fragment), UniProt P23571
Glx (Z), 14-3-3 protein gamma-1, UniProt Q6UFZ3


*/
public class NonstandardProteinCompoundTest  {

	@Test
	public void test() {
		
		
		String[] uniprotIDs = new String[]{"O30642","P59797","Q9N2D9","Q9N2D9","Q6UFZ3"};
		
		for ( String id: uniprotIDs){
			try {
				testUniprot(id);
			} catch (Exception e){
				e.printStackTrace();
				fail(e.getMessage());
				
			}
		}
	}
	
	public void testUniprot(String uniprotID) throws Exception{
		
		ProteinSequence seq = getUniprot(uniprotID);
		
		AminoAcidCompoundSet compoundSet = AminoAcidCompoundSet.getAminoAcidCompoundSet();
		
/*		for (AminoAcidCompound compound : seq) {
			System.out.println(compound.getShortName() + " " + compound.getLongName() + " " + compound.getDescription() + " | " + compoundSet.getEquivalentCompounds(compound) + " " + compound.getMolecularWeight() + " " + compound.getBase());
		} 
		*/
		compoundSet.verifySequence(seq);
		
		
		
		float[] values = Jronn.getDisorderScores(seq);
		
		
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
