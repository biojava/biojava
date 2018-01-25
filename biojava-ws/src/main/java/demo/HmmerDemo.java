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
import org.biojava.nbio.ws.hmmer.HmmerDomain;
import org.biojava.nbio.ws.hmmer.HmmerResult;
import org.biojava.nbio.ws.hmmer.RemoteHmmerScan;

import java.util.SortedSet;

/** 
 * The cookbook recipe for how to request Pfam annotations for a protein sequence using the Hmmer3 service
 *
 * @author Andreas Prlic
 * @since 3.0.3
 */
public class HmmerDemo {

	public static void main(String[] args) throws Exception {


		// first we get a UniProt sequence
		String uniProtID = "P08487";
		ProteinSequence seq = getUniprot(uniProtID);


		// now we submit this sequence to the Hmmer web site
		RemoteHmmerScan hmmer = new RemoteHmmerScan();

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
						domain.getEvalue(), hmmerResult.getDesc()
						));

			}

		}


	}

	/** 
	 * Fetch a protein sequence from the UniProt web site
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
