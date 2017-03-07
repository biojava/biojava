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
package org.biojava.nbio.phylo;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import org.biojava.nbio.core.sequence.MultipleSequenceAlignment;
import org.biojava.nbio.core.sequence.io.FastaWriter;
import org.biojava.nbio.core.sequence.io.template.FastaHeaderFormatInterface;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;
import org.forester.io.parsers.FastaParser;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.msa.Msa;
import org.forester.phylogeny.Phylogeny;

/**
 * This class contains wrapper methods for communication between BioJava and
 * forester (e.g, Data Structure conversion).
 *
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class ForesterWrapper {

	/** Prevent instantiation */
	private ForesterWrapper() {
	}

	/**
	 * Convert a BioJava {@link MultipleSequenceAlignment} to a forester
	 * {@link Msa}. The easiest way to convert them is writting the msa as a
	 * FASTA file and then parsing it with the forester {@link FastaParser}.
	 *
	 * @param msa
	 *            BioJava MultipleSequenceAlignment
	 * @return forester Msa object
	 * @throws IOException
	 *             if the conversion was not possible
	 */
	public static <C extends Sequence<D>, D extends Compound> Msa convert(
			MultipleSequenceAlignment<C, D> msa) throws IOException {

		// Convert the biojava MSA to a FASTA String
		OutputStream os = new ByteArrayOutputStream();
		FastaWriter<C, D> fastaW = new FastaWriter<C, D>(os,
				msa.getAlignedSequences(),
				new FastaHeaderFormatInterface<C, D>() {
					@Override
					public String getHeader(C sequence) {
						return sequence.getAccession().toString();
					};
				});

		fastaW.process();
		String fastaMSA = os.toString();

		// Parse the FASTA file in forester
		return FastaParser.parseMsa(fastaMSA);
	}

	/**
	 * Convert a Phylogenetic tree to its Newick representation, so that it can
	 * be exported to an external application.
	 *
	 * @param phylo
	 *            Phylogeny phylogenetic tree
	 * @param writeDistances
	 *            write the branch lengths if true
	 * @return
	 * @throws IOException
	 */
	public static String getNewickString(Phylogeny phylo,
			boolean writeDistances) throws IOException {

		PhylogenyWriter w = new PhylogenyWriter();
		StringBuffer newickString = w.toNewHampshire(phylo, writeDistances);
		return newickString.toString();
	}

	/**
	 * Helper function to clone a forester symmetrical DistanceMatrix.
	 *
	 * @param distM
	 *            forester symmetrical DistanceMatrix
	 * @return identical copy of the forester symmetrical DistanceMatrix
	 */
	public static BasicSymmetricalDistanceMatrix cloneDM(
			BasicSymmetricalDistanceMatrix distM) {

		int n = distM.getSize();
		BasicSymmetricalDistanceMatrix cloneDM =
				new BasicSymmetricalDistanceMatrix(n);

		for (int i = 0; i < n; i++) {
			cloneDM.setIdentifier(i, distM.getIdentifier(i));
			for (int j = i + 1; j < n; j++) {
				cloneDM.setValue(i, j, distM.getValue(i, j));
			}
		}
		return cloneDM;
	}

}
