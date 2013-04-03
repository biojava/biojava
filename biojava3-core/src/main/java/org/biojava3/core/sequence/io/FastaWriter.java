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
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.Collection;
import java.util.LinkedHashMap;

import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.io.template.FastaHeaderFormatInterface;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

/**
 * The FastaWriter writes a collection of sequences to an outputStream. FastaWriterHelper should be
 * used to write out sequences. Each sequence loaded from a fasta file retains the original Fasta header
 * and that is used when writing to the stream. This behavior can be overwritten by implementing
 * a custom FastaHeaderFormatInterface.
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class FastaWriter<S extends Sequence<?>, C extends Compound> {

    OutputStream os;
    Collection<S> sequences;
    FastaHeaderFormatInterface<S, C> headerFormat;
    private int lineLength = 60;
    byte[] lineSep = System.getProperty("line.separator").getBytes();
/**
 * Use default line length of 60
 * @param os
 * @param sequences
 * @param headerFormat
 */
    public FastaWriter(OutputStream os, Collection<S> sequences, FastaHeaderFormatInterface<S, C> headerFormat) {

        this.os = os;
        this.sequences = sequences;
        this.headerFormat = headerFormat;
    }

/**
 * Set custom lineLength
 * @param os
 * @param sequences
 * @param headerFormat
 * @param lineLength
 */

    public FastaWriter(OutputStream os, Collection<S> sequences, FastaHeaderFormatInterface<S, C> headerFormat, int lineLength) {
        this.os = os;
        this.sequences = sequences;
        this.headerFormat = headerFormat;
        this.lineLength = lineLength;
    }

    /**
     * Allow an override of operating system line separator for programs that needs a specific CRLF or CR or LF option
     * @param lineSeparator
     */
    public void setLineSeparator(String lineSeparator){
        lineSep = lineSeparator.getBytes();
    }

    public void process() throws Exception {
       // boolean closeit = false;
       
        

        for (S sequence : sequences) {
            String header = headerFormat.getHeader(sequence);
            os.write('>');
            os.write(header.getBytes());
            os.write(lineSep);

            int compoundCount = 0;
            String seq = "";

            seq = sequence.getSequenceAsString();

            for (int i = 0; i < seq.length(); i++) {
                os.write(seq.charAt(i));
                compoundCount++;
                if (compoundCount == lineLength) {
                    os.write(lineSep);
                    compoundCount = 0;
                }

            }


            //If we had sequence which was a reciprocal of line length
            //then don't write the line terminator as this has already written
            //it
            if ((sequence.getLength() % getLineLength()) != 0) {
                os.write(lineSep);
            }
        }
        
    }

    public static void main(String[] args) {
        try {
            FileInputStream is = new FileInputStream("/Users/Scooter/scripps/dyadic/c1-454Scaffolds.faa");


            FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(is, new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(), new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
            LinkedHashMap<String, ProteinSequence> proteinSequences = fastaReader.process();
            is.close();


          //  System.out.println(proteinSequences);

            FileOutputStream fileOutputStream = new FileOutputStream("/Users/Scooter/scripps/dyadic/c1-454Scaffolds_temp.faa");
           
            BufferedOutputStream bo = new BufferedOutputStream(fileOutputStream);
            long start = System.currentTimeMillis();
            FastaWriter<ProteinSequence, AminoAcidCompound> fastaWriter = new FastaWriter<ProteinSequence, AminoAcidCompound>(bo, proteinSequences.values(), new GenericFastaHeaderFormat<ProteinSequence, AminoAcidCompound>());
            fastaWriter.process();
            bo.close();
            long end = System.currentTimeMillis();
            System.out.println("Took " + (end - start) + " seconds");
         
            fileOutputStream.close();


        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * @return the lineLength
     */
    public int getLineLength() {
        return lineLength;
    }

    /**
     * @param lineLength the lineLength to set
     */
    public void setLineLength(int lineLength) {
        this.lineLength = lineLength;
    }
}
