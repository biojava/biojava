package demo;

import org.biojava.nbio.core.sequence.MultipleSequenceAlignment;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.phylo.TreeConstructionAlgorithm;
import org.biojava.nbio.phylo.ProgressListenerStub;
import org.biojava.nbio.phylo.TreeConstructor;
import org.biojava.nbio.phylo.TreeType;

/**
 * This demo contains the CookBook examples to create a phylogenetic tree from a
 * given multiple sequence alignment (MSA).
 * 
 * @author Aleix Lafita
 *
 */
public class DemoTreeConstructor {

	public static void main(String[] args) throws Exception {

		MultipleSequenceAlignment<ProteinSequence, AminoAcidCompound> msa = new 
				MultipleSequenceAlignment<ProteinSequence, AminoAcidCompound>();

		ProteinSequence pSeq1 = new ProteinSequence(
				"RER-RDGGGNSRKYDDRRSPRDGE---IDYDERTVSHYQRQFQDERISDGM"
				+ "LNTLKQSLKGLDCQPIHLKDSKANRSIMIDEIHTGTADSVTFEQKLPDGEMKL");
		ProteinSequence pSeq2 = new ProteinSequence(
				"RDRHRD---DRHRYDEDRDHRRDQRNVSDYDSEELRKFEEDYKSDRLGQYV"
				+ "FSDLNSAVKGLVVQPIHL-NKEVNRTVIIDSICKESAEKVRFEFGKGEDAREI");
		ProteinSequence pSeq3 = new ProteinSequence(
				"RPTH---GGLSLNIDVSTTMILEPGPVIEF-----LKANQSVETPRQIDWI"
				+ "-----KAAKML--KHMRVKATHRNMEFKIIGLSSKPCNQQLFSMKIKDGEREV");

		msa.addAlignedSequence(pSeq1);
		msa.addAlignedSequence(pSeq2);
		msa.addAlignedSequence(pSeq3);

		TreeConstructor<ProteinSequence, AminoAcidCompound> treeConstructor = 
				new TreeConstructor<ProteinSequence, AminoAcidCompound>(
				msa, TreeType.NJ, TreeConstructionAlgorithm.PID, 
				new ProgressListenerStub());
		
		treeConstructor.process();
		String newick = treeConstructor.getNewickString(true, true);

		System.out.println(newick);
	}
}
