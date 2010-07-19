/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence.io;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collection;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.GeneSequence;
import org.biojava3.core.sequence.ProteinSequence;

import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.template.FastaHeaderFormatInterface;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class FastaWriterHelper {



    /**
     * Write collection of protein sequences to a file
     *
     * @param file
     * @param proteinSequences
     * @throws Exception
     */
    public static void writeProteinSequence(File file,
            Collection<ProteinSequence> proteinSequences) throws Exception {
        FileOutputStream outputStream = new FileOutputStream(file);
        writeProteinSequence(outputStream, proteinSequences);
        outputStream.close();
    }

    /**
     * Write collection of protein sequences to a stream
     * @param outputStream
     * @param proteinSequences
     * @throws Exception
     */

    public static void writeProteinSequence(OutputStream outputStream,
            Collection<ProteinSequence> proteinSequences) throws Exception {
        FastaWriter<ProteinSequence, AminoAcidCompound> fastaWriter = new FastaWriter<ProteinSequence, AminoAcidCompound>(
                outputStream, proteinSequences,
                new GenericFastaHeaderFormat<ProteinSequence, AminoAcidCompound>());
        fastaWriter.process();

    }

        /**
     * Write a collection of GeneSequences to a file where if the gene is negative strand it will flip and complement the sequence
     * @param file
     * @param geneSequences
     * @throws Exception
     */

    public static void writeGeneSequence(File file, Collection<GeneSequence> geneSequences,boolean showExonUppercase) throws Exception {
        FileOutputStream outputStream = new FileOutputStream(file);
        writeGeneSequence(outputStream, geneSequences,showExonUppercase);
        outputStream.close();
    }

    /**
     * Write a collection of GeneSequences to a file where if the gene is negative strand it will flip and complement the sequence
     * @param outputStream
     * @param dnaSequences
     * @throws Exception
     */

    public static void writeGeneSequence(OutputStream outputStream, Collection<GeneSequence> geneSequences,boolean showExonUppercase) throws Exception {
        FastaGeneWriter fastaWriter = new FastaGeneWriter(
                outputStream, geneSequences,
                new GenericFastaHeaderFormat<GeneSequence, NucleotideCompound>(),showExonUppercase);
        fastaWriter.process();

    }


    /**
     * Write a collection of NucleotideSequences to a file
     * @param file
     * @param dnaSequences
     * @throws Exception
     */

    public static void writeNucleotideSequence(File file, Collection<DNASequence> dnaSequences) throws Exception {
        FileOutputStream outputStream = new FileOutputStream(file);
        writeNucleotideSequence(outputStream, dnaSequences);
        outputStream.close();
    }

    /**
     * Write a collection of NucleotideSequences to a file
     * @param outputStream
     * @param dnaSequences
     * @throws Exception
     */

    public static void writeNucleotideSequence(OutputStream outputStream, Collection<DNASequence> dnaSequences) throws Exception {
        FastaWriter<DNASequence, NucleotideCompound> fastaWriter = new FastaWriter<DNASequence, NucleotideCompound>(
                outputStream, dnaSequences,
                new GenericFastaHeaderFormat<DNASequence, NucleotideCompound>());
        fastaWriter.process();

    }

    public static void writeSequence(File file, Sequence<?> sequence) throws Exception {
        writeSequences(new FileOutputStream(file), singleSeqToCollection(sequence));
    }

    public static void writeSequence(OutputStream outputStream, Sequence<?> sequence) throws Exception {
        writeSequences(outputStream, singleSeqToCollection(sequence));
    }

    private static Collection<Sequence<?>> singleSeqToCollection(Sequence<?> sequence) {
        Collection<Sequence<?>> sequences = new ArrayList<Sequence<?>>();
        sequences.add(sequence);
        return sequences;
    }

    /**
     * Method which will write your given Sequences to the specified
     * {@link OutputStream}. This is a very generic method which writes just the
     * AccessionID of the Sequence as the FASTA header.
     *
     * @param outputStream Stream to write to; can be System.out
     * @param sequences The sequences to write out
     * @throws Exception Thrown normally thanks to IO problems
     */
    public static void writeSequences(OutputStream outputStream,
            Collection<Sequence<?>> sequences) throws Exception {

        FastaHeaderFormatInterface<Sequence<?>, Compound> fhfi =
                new FastaHeaderFormatInterface<Sequence<?>, Compound>() {

                    public String getHeader(Sequence<?> sequence) {
                        return sequence.getAccession().toString();
                    }

                    ;
                };

        FastaWriter<Sequence<?>, Compound> fastaWriter =
                new FastaWriter<Sequence<?>, Compound>(outputStream,
                sequences, fhfi);

        fastaWriter.process();
    }
}
