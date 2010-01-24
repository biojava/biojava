/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence.io;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.Collection;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.io.template.FastaHeaderFormatInterface;
import org.biojava3.core.sequence.template.AbstractSequence;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class FastaWriter<S extends AbstractSequence> {

    OutputStream os;
    Collection<S> sequences;
    FastaHeaderFormatInterface headerFormat;
    private int lineLength = 60;

    public FastaWriter(OutputStream os, Collection<S> sequences, FastaHeaderFormatInterface headerFormat) {
        this.os = os;
        this.sequences = sequences;
        this.headerFormat = headerFormat;
    }

    public void process() throws Exception {
        for (AbstractSequence sequence : sequences) {
            String header = headerFormat.getHeader(sequence);
            os.write('>');
            os.write(header.getBytes());
            os.write('\n');
            byte[] data = sequence.getSequenceAsString().getBytes();

            int bytesToWrite = 0;
            for(int i = 0; i < data.length; i = i + getLineLength()){
                if(i + getLineLength() < data.length){
                  bytesToWrite = getLineLength();
                }else{
                  bytesToWrite = data.length - i;
                }
              
                os.write(data, i,bytesToWrite);
                os.write('\n');
            }
            
        }
    }

    public static void main(String[] args) {
        try {
            FileInputStream is = new FileInputStream("/Users/Scooter/mutualinformation/project/nuclear_receptor/PF00104_small.fasta");
            

            FastaReader<ProteinSequence> fastaReader = new FastaReader<ProteinSequence>(is, new GenericFastaHeaderParser(), new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
            Collection<ProteinSequence> proteinSequences = fastaReader.process();
            is.close();
            

            System.out.println(proteinSequences);
            
            FileOutputStream fileOutputStream = new FileOutputStream("/Users/Scooter/mutualinformation/project/nuclear_receptor/PF00104_small_test.fasta");
            
            FastaWriter fastaWriter = new FastaWriter(fileOutputStream,proteinSequences,new GenericFastaHeaderFormat());
            fastaWriter.process();
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
