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
import org.biojava.nbio.phylo.DistanceMatrixCalculator;
import org.biojava.nbio.phylo.DistanceTreeEvaluator;
import org.biojava.nbio.phylo.ForesterWrapper;
import org.biojava.nbio.phylo.TreeConstructor;
import org.biojava.nbio.phylo.TreeConstructorType;
import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;
import org.forester.evoinference.matrix.distance.DistanceMatrix;
import org.forester.phylogeny.Phylogeny;

/**
 * This demo contains the CookBook example to create a phylogenetic tree from a
 * given multiple sequence alignment (MSA).
 * 
 * @author Scooter Willis
 * @author Aleix Lafita
 *
 */
public class DemoDistanceTree {

	public static void main(String[] args) throws Exception {

		// 0. This is just to load an example MSA from a FASTA file
		InputStream inStream = TreeConstructor.class
				.getResourceAsStream("/PF00104_small.fasta");

		FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(
				inStream,
				new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(),
				new ProteinSequenceCreator(AminoAcidCompoundSet
						.getAminoAcidCompoundSet()));

		LinkedHashMap<String, ProteinSequence> proteinSequences = fastaReader
				.process();

		inStream.close();

		MultipleSequenceAlignment<ProteinSequence, AminoAcidCompound> msa = new MultipleSequenceAlignment<ProteinSequence, AminoAcidCompound>();

		for (ProteinSequence proteinSequence : proteinSequences.values()) {
			msa.addAlignedSequence(proteinSequence);
		}

		long readT = System.currentTimeMillis();

		// 1. Calculate the evolutionary distance matrix
		DistanceMatrix DM = DistanceMatrixCalculator.kimuraDistance(msa);

		// 2. Construct a distance tree using the NJ algorithm
		Phylogeny phylo = TreeConstructor.distanceTree(
				(BasicSymmetricalDistanceMatrix) DM, TreeConstructorType.NJ);

		long treeT = System.currentTimeMillis();
		String newick = ForesterWrapper.getNewickString(phylo, true, true);
		System.out.println(newick);
		System.out.println("Tree Construction: " + (treeT - readT) + " ms.");

		// 3. Evaluate the goodness of fit of the tree
		double cv = DistanceTreeEvaluator.evaluate(phylo, DM);
		System.out.println("CV of the tree: " + (int) (cv*100) + " %");

	}
}
