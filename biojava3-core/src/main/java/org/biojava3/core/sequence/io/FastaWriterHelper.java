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
 * Created on 01-21-2010
 */
package org.biojava3.core.sequence.io;

import java.io.BufferedOutputStream;
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
 * The class that should be used to write out fasta file of a sequence collection
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
        BufferedOutputStream bo = new BufferedOutputStream(outputStream);
        writeProteinSequence(bo, proteinSequences);
        bo.close();
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
        BufferedOutputStream bo = new BufferedOutputStream(outputStream);
        writeGeneSequence(bo, geneSequences,showExonUppercase);
        bo.close();
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
        BufferedOutputStream bo = new BufferedOutputStream(outputStream);
        writeNucleotideSequence(bo, dnaSequences);
        bo.close();
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

    /**
     * Write a sequence to a file
     * @param file
     * @param sequence
     * @throws Exception
     */
    public static void writeSequence(File file, Sequence<?> sequence) throws Exception {
        FileOutputStream outputStream = new FileOutputStream(file);
        BufferedOutputStream bo = new BufferedOutputStream(outputStream);
        writeSequences(bo, singleSeqToCollection(sequence));
        bo.close();
        outputStream.close();
    }

    /**
     * Write a sequence to OutputStream
     * @param outputStream
     * @param sequence
     * @throws Exception
     */
    public static void writeSequence(OutputStream outputStream, Sequence<?> sequence) throws Exception {
        writeSequences(outputStream, singleSeqToCollection(sequence));
    }

    /**
     * 
     * @param sequence
     * @return
     */

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
