package demo;

import java.util.SortedSet;


import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.loader.UniprotProxySequenceReader;
import org.biojava3.ws.hmmer.HmmerDomain;
import org.biojava3.ws.hmmer.HmmerResult;
import org.biojava3.ws.hmmer.HmmerScan;
import org.biojava3.ws.hmmer.RemoteHmmerScan;

/** The cookbook recipe for how to request Pfam annotations for a protein sequence using the Hmmer3 service
 * 
 * @author Andreas Prlic
 * @since 3.0.3
 */
public class HmmerDemo {

	public static void main(String[] args){

		try {
			// first we get a UniProt sequence
			String uniProtID = "P26663";
			ProteinSequence seq = getUniprot(uniProtID);

			
			// now we submit this sequence to the Hmmer web site
			HmmerScan hmmer = new RemoteHmmerScan();

			SortedSet<HmmerResult> results = hmmer.scan(seq);

			
			
			// and now let's print out the obtained annotations
			
			System.out.println(String.format("#\t%15s\t%10s\t%s\t%s\t%8s\t%s",
					"Domain","ACC", "Start","End","eValue","Description"));
			
			int counter = 0;
			for (HmmerResult hmmerResult : results) {
				//System.out.println(hmmerResult);

				for ( HmmerDomain domain : hmmerResult.getDomains()) {
					counter++;
					System.out.println(String.format("%d\t%15s\t%10s\t%5d\t%5d\t%.2e\t%s",
							counter,
							hmmerResult.getName(), domain.getHmmAcc(), 
							domain.getSqFrom(),domain.getSqTo(),
							hmmerResult.getEvalue(), hmmerResult.getDesc()
							));

				}

			}

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
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
