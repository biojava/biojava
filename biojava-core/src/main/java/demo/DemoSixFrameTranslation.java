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

import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.AmbiguityRNACompoundSet;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.io.DNASequenceCreator;
import org.biojava.nbio.core.sequence.io.FastaReader;
import org.biojava.nbio.core.sequence.io.GenericFastaHeaderParser;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.core.sequence.transcription.Frame;
import org.biojava.nbio.core.sequence.transcription.TranscriptionEngine;

import java.io.ByteArrayInputStream;
import java.io.InputStream;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Created by andreas on 8/10/15.
 */
public class DemoSixFrameTranslation {

	public static void main(String[] args){
		String dnaFastaS = ">gb:GQ903697|Organism:Arenavirus H0030026 H0030026|Segment:S|Host:Rat\n" +
				"CGCACAGAGGATCCTAGGCGTTACTGACTTGCGCTAATAACAGATACTGTTTCATATTTAGATAAAGACC\n" +
				"CAGCCAACTGATTGGTCAGCATGGGACAACTTGTGTCCCTCTTCAGTGAAATTCCATCAATCATACACGA\n" +
				"AGCTCTCAATGTTGCTCTCGTAGCTGTTAGCATCATTGCAATATTGAAAGGGGTTGTGAATGTTTGGAAG\n" +
				"AGTGGAGTTTTGCAGCTTTTGGCCTTCTTGCTCCTGGCGGGAAGATCCTGCTCAGTCATAATTGGTCATC\n" +
				"ATCTCGAACTGCAGCATGTGATCTTCAATGGGTCATCAATCACACCCTTTTTACCAGTTACATGTAAGAT\n" +
				"CAATGATACCTACTTCCTACTAAGAGGCCCCTATGAAGCTGATTGGGCAGTTGAATTGAGTGTAACTGAA\n" +
				"ACCACAGTCTTGGTTGATCTTGAAGGTGGCAGCTCAATGAAGCTGAAAGCCGGAAACATCTCAGGTTGTC\n" +
				"TTGGAGACAACCCCCATCTGAGATCAGTGGTCTTCACATTGAATTGGTTGCTAACAGGATTAGATCATGT\n" +
				"TATTGATTCTGACCCGAAAATTCTCTGTGATCTTAAAGACAGTGGGCACTTTCGTCTCCAGATGAACTTA\n" +
				"ACAGAAAAGCACTATTGTGACAAGTTTCACATCAAAATGGGCAAGGTCTTTGGCGTATTCAAAGATCCGT\n" +
				"GCATGGCTGGTGGTAAAATGTTTGCCATACTAAAAAATACCTCTTGGTCGAACCAGTGCCAAGGAAACCA\n" +
				"TGTCAGCACCATTCATCTTGTCCTTCAGAGTAATTTCAAACAGGTCCTCAGTAGCAGGAAACTGTTGAAC\n" +
				"TTTTTCAGCTGGTCATTGTCTGATGCCACAGGGGCTGATATGCCTGGTGGTTTTTGTCTGGAAAAATGGA\n" +
				"TGTTGATTTCAAGTGAACTGAAATGCTTTGGAAACACAGCTGTGGCAAAGTGCAACTTAAATCATGACTC\n" +
				"AGAGTTCTGTGACATGCTTAGGCTTTTTGATTTCAACAAAAAGGCAATAGTCACTCTTCAGAACAAAACA\n" +
				"AAGCATCGGCTGGACACAGTAATTACTGCTATCAATTCATTGATCTCTGATAATATTCTTATGAAGAACA\n" +
				"GGATTAAAGAATTGATAGATGTTCCTTACTGTAATTACACCAAATTTTGGTATGTCAATCACACAGGTCT\n" +
				"AAATCTGCACACCCTTCCAAGATGTTGGCTTGTTAAAAATGGTAGCTACTTGAATGTGTCTGACTTCAGG\n" +
				"AATGAGTGGATATTGGAGAGTGATCATCTTGTTTCGGAGATCCTTTCAAAGGAGTATGAGGAAAGGCAAA\n" +
				"ATCGTACACCACTCTCACTGGTTGACATCTGTTTCTGGAGTACATTGTTTTACACAGCATCAATTTTCCT\n" +
				"ACACCTCTTGAGAATTCCAACCCACAGACACATTGTTGGTGAGGGCTGCCCGAAGCCTCATAGGCTAAAC\n" +
				"AGGCACTCAATATGTGCTTGTGGCCTTTTCAAACAAGAAGGCAGACCCTTGAGATGGGTAAGAAAGGTGT\n" +
				"GAACAATGGTTGCTTGGTGGCCTCCATTGCTGCACCCCCCTAGGGGGGTGCAGCAATGGAGGTTCTCGYT\n" +
				"GAGCCTAGAGAACAACTGTTGAATCGGGTTCTCTAAAGAGAACATCGATTGGTAGTACCCTTTTTGGTTT\n" +
				"TTCATTGGTCACTGACCCTGAAAGCACAGCACTGAACATCAAACAGTCCAAAAGTGCACAGTGTGCATTT\n" +
				"GTTGTGGCTGGTGCTGATCCTTTCTTCTTACTTTTAATGACTATTCCCTTATGTCTGTCACACAGATGTT\n" +
				"CAAATCTCTTCCAAACAAGATCTTCAAAGAGCCGTGACTGTTCTGCGGTCAGTTTGACATCAACAATCTT\n" +
				"CAAATCCTGTCTTCCATGCATATCAAAGAGCCTCCTAATATCATCAGCACCTTGCGCAGTGAAAACCATG\n" +
				"GATTTAGGCAGACTCCTTATTATGCTTGTGATGAGGCCAGGTCGTGCATGTTCAACATCCTTCAGCAATA\n" +
				"TCCCATGACAATATTTACTTTGGTCCTTAAAAGATTTTATGTCATTGGGTTTTCTGTAGCAGTGGATGAA\n" +
				"TTTTTGTGATTCAGGCTGGTAAATTGCAAACTCAACAGGGTCATGTGGCGGGCCTTCAATGTCAATCCAT\n" +
				"GTTGTGTCACTGACCATCAACGACTCTACACTTCTCTTCACCTGAGCCTCCACCTCAGGCTTGAGCGTGG\n" +
				"ACAAGAGTGGGGCACCACCGTTCCGGATGGGGACTGGTGTTTTGCTTGGTAAACTCTCAAATTCCACAAC\n" +
				"TGTATTGTCCCATGCTCTCCCTTTGATCTGTGATCTTGATGAAATGTAAGGCCAGCCCTCACCAGAGAGA\n" +
				"CACACCTTATAAAGTATGTTTTCATAAGGATTCCTCTGTCCTGGTATGGCACTGATGAACATGTTTTCCC\n" +
				"TCTTTTTGATCTCCAAGAGGGTTTTTATAATGGTTGTGAATGTGGACTCCTCAATCTTTATTGTTTCCAG\n" +
				"CATGTTGCCACCATCAATCAGGCAAGCACCGGCTTTCACAGCAGCTGATAAACTAAGGTTGTAGCCTGAT\n" +
				"ATGTTAATTTGAGAATCCTCCTGAGTGATTACCTTTAGAGAAGGATGCTTCTCCATCAAAGCATCTAAGT\n" +
				"CACTTAAATTAGGGTATTTTGCTGTGTATAGCAACCCCAGATCTGTGAGGGCCTGAACCACATCATTTAG\n" +
				"AGTTTCCCCTCCCTGTTCAGTCATACAGGAAATTGTGAGTGCTGGCATCGATCCAAATTGGTTGATCATA\n" +
				"AGTGATGAGTCTTTAACGTCCCAGACTTTGACCACCCCTCCAGTTCTAGCCAACCCAGGTCTCTGAATAC\n" +
				"CAACAAGTTGCAGAATTTCGGACCTCCTGGTGAGCTGTGTTGTAGAGAGGTTCCCTAGATACTGGCCACC\n" +
				"TGTGGCTGTCAACCTCTCTGTTCTTTGAACTTTTTGCCTTAATTTGTCCAAGTCACTGGAGAGTTCCATT\n" +
				"AGCTCTTCCTTTGACAATGATCCTATCTTAAGGAACATGTTCTTTTGGGTTGACTTCATGACCATCAATG\n" +
				"AGTCAACTTCCTTATTCAAGTCCCTCAAACTAACAAGATCACTGTCATCTCTTTTAGACCTCCTCATCAT\n" +
				"GCGTTGCACACTTGCAACCTTTGAAAAATCTAAGCCGGACAGAAGAGCCCTCGCGTCAGTTAGGACATCT\n" +
				"GCCTTAACAGCAGTTGTCCAGTTCGAGAGTCCTCTCCTGAGAGACTGTGTCCATCTGAATGATGGGATTG\n" +
				"GTTGTTCGCTCATAGTGATGAAATTGCGCAGAGTTATCCAAAAGCCTAGGATCCTCTGTGCG";


		try {

			// parse the raw sequence from the string
			InputStream stream = new ByteArrayInputStream(dnaFastaS.getBytes());

			// define the Ambiguity Compound Sets
			AmbiguityDNACompoundSet ambiguityDNACompoundSet = AmbiguityDNACompoundSet.getDNACompoundSet();
			CompoundSet<NucleotideCompound> nucleotideCompoundSet = AmbiguityRNACompoundSet.getRNACompoundSet();

			FastaReader<DNASequence, NucleotideCompound> proxy =
					new FastaReader<DNASequence, NucleotideCompound>(
							stream,
							new GenericFastaHeaderParser<DNASequence, NucleotideCompound>(),
							new DNASequenceCreator(ambiguityDNACompoundSet));

			// has only one entry in this example, but could be easily extended to parse a FASTA file with multiple sequences
			LinkedHashMap<String, DNASequence> dnaSequences = proxy.process();

			// Initialize the Transcription Engine
			TranscriptionEngine engine = new
					TranscriptionEngine.Builder().dnaCompounds(ambiguityDNACompoundSet).rnaCompounds(nucleotideCompoundSet).build();

			Frame[] sixFrames = Frame.getAllFrames();



			for (DNASequence dna : dnaSequences.values()) {

				Map<Frame, Sequence<AminoAcidCompound>> results = engine.multipleFrameTranslation(dna, sixFrames);

				for (Frame frame : sixFrames){
					System.out.println("Translated Frame:" + frame +" : " + results.get(frame));
					//System.out.println(dna.getRNASequence(frame).getProteinSequence(engine));

					ProteinSequence ps = new ProteinSequence(results.get(frame).getSequenceAsString());
					System.out.println(ps);
					try {

					} catch (Exception e){
						System.err.println(e.getMessage() + " when trying to translate frame " + frame);
					}
				}

			}
		} catch (Exception e){
			e.printStackTrace();
		}


	}



}
