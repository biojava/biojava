/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence.io;

import java.io.File;
import java.io.FileInputStream;
import org.biojava3.core.sequence.io.template.FastaHeaderParserInterface;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.template.AbstractSequence;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class FastaReader<S extends AbstractSequence> {

    SequenceCreatorInterface sequenceCreator;
    FastaHeaderParserInterface headerParser;
    BufferedReaderBytesRead br;
    InputStreamReader isr;
    FileInputStream fi = null;

    /**
     *
     * @param br
     * @param headerParser
     * @param sequenceCreator
     */
    public FastaReader(InputStream is, FastaHeaderParserInterface headerParser, SequenceCreatorInterface sequenceCreator) {
        this.headerParser = headerParser;
        isr = new InputStreamReader(is);
        this.br = new BufferedReaderBytesRead(isr);
        this.sequenceCreator = sequenceCreator;
    }

    public FastaReader(File file, FastaHeaderParserInterface headerParser, SequenceCreatorInterface sequenceCreator) throws Exception {
        fi = new FileInputStream(file);
        isr = new InputStreamReader(fi);
        this.br = new BufferedReaderBytesRead(isr);
        this.sequenceCreator = sequenceCreator;
    }

    /**
     * The parsing is done in this method
     * @return
     * @throws Exception
     */
    public List<S> process() throws Exception {
        ArrayList<S> sequences = new ArrayList<S>();


        String line = "";
        String header = "";
        StringBuffer sb = new StringBuffer();
        int maxSequenceLength = -1;
        long fileIndex = 0;
        long sequenceIndex = 0;
        boolean keepGoing = true;
        do {
            line = line.trim(); // nice to have but probably not needed
            if (line.length() != 0) {
                if (line.startsWith(">")) {
                    if (sb.length() > 0) {
                        S sequence = (S) sequenceCreator.getSequence(sb.toString(), sequenceIndex);
                        headerParser.parseHeader(header, sequence);
                        sequences.add(sequence);
                        if (maxSequenceLength < sb.length()) {
                            maxSequenceLength = sb.length();
                        }
                        sb = new StringBuffer(maxSequenceLength);
                    }
                    header = line.substring(1);
                } else if (line.startsWith(";")) {
                } else {
                    //mark the start of the sequence with the fileIndex before the line was read
                    if(sb.length() == 0){
                        sequenceIndex = fileIndex;
                    }
                    sb.append(line);
                }
            }
            fileIndex = br.getBytesRead();
            line = br.readLine();
            if (line == null) {
                S sequence = (S) sequenceCreator.getSequence(sb.toString(), fileIndex);
                headerParser.parseHeader(header, sequence);
                sequences.add(sequence);
                keepGoing = false;
            }
        } while (keepGoing);
        br.close();
        isr.close();
        //If stream was created from File object then we need to close it
        if (fi != null) {
            fi.close();
        }
        return sequences;
    }

    public static void main(String[] args) {
        try {
            FileInputStream is = new FileInputStream("/Users/Scooter/mutualinformation/project/nuclear_receptor/PF00104_small.fasta");


            FastaReader<ProteinSequence> fastaReader = new FastaReader<ProteinSequence>(is, new GenericFastaHeaderParser(), new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
            Collection<ProteinSequence> proteinSequences = fastaReader.process();
            is.close();


            System.out.println(proteinSequences);

            File file = new File("/Users/Scooter/mutualinformation/project/nuclear_receptor/PF00104_small.fasta");
            FastaReader<ProteinSequence> fastaProxyReader = new FastaReader<ProteinSequence>(file, new GenericFastaHeaderParser(), new FileProxyProteinSequenceCreator(file, AminoAcidCompoundSet.getAminoAcidCompoundSet()));
            Collection<ProteinSequence> proteinProxySequences = fastaProxyReader.process();

            System.out.println(proteinProxySequences);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
