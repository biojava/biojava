package demo;

import java.io.InputStream;
import java.util.LinkedHashMap;

import org.biojava.nbio.core.sequence.MultipleSequenceAlignment;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.io.FastaReader;
import org.biojava.nbio.core.sequence.io.GenericFastaHeaderParser;
import org.biojava.nbio.core.sequence.io.ProteinSequenceCreator;
import org.biojava.nbio.phylo.TreeConstructionAlgorithm;
import org.biojava.nbio.phylo.ProgressListenerStub;
import org.biojava.nbio.phylo.TreeConstructor;
import org.biojava.nbio.phylo.TreeType;

/**
 * This demo contains the CookBook example to create a phylogenetic tree from a
 * given multiple sequence alignment (MSA).
 * 
 * @author Scooter Willis
 * @author Aleix Lafita
 *
 */
public class DemoTreeConstructor {

	public static void main(String[] args) throws Exception {

		InputStream inStream = TreeConstructor.class
				.getResourceAsStream("/PF00104_small.fasta");

		FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = 
				new FastaReader<ProteinSequence, AminoAcidCompound>(
				inStream, new GenericFastaHeaderParser<ProteinSequence, 
				AminoAcidCompound>(), new ProteinSequenceCreator(
						AminoAcidCompoundSet.getAminoAcidCompoundSet()));
		
		LinkedHashMap<String, ProteinSequence> proteinSequences = 
				fastaReader.process();
		
		inStream.close();

		MultipleSequenceAlignment<ProteinSequence, AminoAcidCompound> msa = 
				new MultipleSequenceAlignment<ProteinSequence, 
				AminoAcidCompound>();
		
		for (ProteinSequence proteinSequence : proteinSequences.values()) {
			msa.addAlignedSequence(proteinSequence);
		}

		long readTime = System.currentTimeMillis();
		TreeConstructor<ProteinSequence, AminoAcidCompound> treeConstructor = 
				new TreeConstructor<ProteinSequence, AminoAcidCompound>(
				msa, TreeType.NJ, TreeConstructionAlgorithm.PID,
				new ProgressListenerStub());
		
		treeConstructor.process();
		long treeTime = System.currentTimeMillis();
		String newick = treeConstructor.getNewickString(true, true);

		System.out.println(newick);
		System.out.println("Tree time {" + (treeTime - readTime) + "}");
		
	}
}
