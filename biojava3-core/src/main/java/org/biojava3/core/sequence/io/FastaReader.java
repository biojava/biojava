/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence.io;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.LinkedHashMap;

import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.io.template.FastaHeaderParserInterface;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class FastaReader<S extends Sequence<?>, C extends Compound> {

    SequenceCreatorInterface<C> sequenceCreator;
    FastaHeaderParserInterface<S,C> headerParser;
    BufferedReaderBytesRead br;
    InputStreamReader isr;
    FileInputStream fi = null;

    /**
     * If you are going to use FileProxyProteinSequenceCreator then do not use this constructor because we need details about
     * local file offsets for quick reads. InputStreams don't give you the name of the stream to access quickly via file seek. A seek in
     * an inputstream is forced to read all the data so you don't gain anything.
     * @param br
     * @param headerParser
     * @param sequenceCreator
     */
    public FastaReader(InputStream is, FastaHeaderParserInterface<S,C> headerParser, SequenceCreatorInterface<C> sequenceCreator) {
        this.headerParser = headerParser;
        isr = new InputStreamReader(is);
        this.br = new BufferedReaderBytesRead(isr);
        this.sequenceCreator = sequenceCreator;
    }

    /**
     * If you are going to use the FileProxyProteinSequenceCreator then you need to use this constructor because we need details about
     * the location of the file.
     * @param file
     * @param headerParser
     * @param sequenceCreator
     * @throws Exception
     */
    public FastaReader(File file, FastaHeaderParserInterface<S,C> headerParser, SequenceCreatorInterface<C> sequenceCreator) throws Exception {
        this.headerParser = headerParser;
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
    @SuppressWarnings("unchecked")
    public LinkedHashMap<String,S> process() throws Exception {
        LinkedHashMap<String,S> sequences = new LinkedHashMap<String,S>();


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
                    //    System.out.println("Sequence index=" + sequenceIndex);
                        S sequence = (S)sequenceCreator.getSequence(sb.toString(), sequenceIndex);
                        headerParser.parseHeader(header, sequence);
                        sequences.put(sequence.getAccession().getID(),sequence);
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
            //    System.out.println("Sequence index=" + sequenceIndex + " " + fileIndex );
                S sequence = (S)sequenceCreator.getSequence(sb.toString(), sequenceIndex);
                headerParser.parseHeader(header, sequence);
                sequences.put(sequence.getAccession().getID(),sequence);
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
            String inputFile = "src/test/resources/PF00104_small.fasta";
            FileInputStream is = new FileInputStream(inputFile);

            FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(is, new GenericFastaHeaderParser<ProteinSequence,AminoAcidCompound>(), new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
            LinkedHashMap<String,ProteinSequence> proteinSequences = fastaReader.process();
            is.close();


            System.out.println(proteinSequences);

            File file = new File(inputFile);
            FastaReader<ProteinSequence,AminoAcidCompound> fastaProxyReader = new FastaReader<ProteinSequence,AminoAcidCompound>(file, new GenericFastaHeaderParser<ProteinSequence,AminoAcidCompound>(), new FileProxyProteinSequenceCreator(file, AminoAcidCompoundSet.getAminoAcidCompoundSet()));
            LinkedHashMap<String,ProteinSequence> proteinProxySequences = fastaProxyReader.process();

            for(String key : proteinProxySequences.keySet()){
                ProteinSequence proteinSequence = proteinProxySequences.get(key);
                System.out.println(key);
                if(key.equals("Q98SJ1_CHICK/15-61")){
                    int dummy = 1;
                }
                System.out.println(proteinSequence.toString());

            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
