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

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.LinkedHashMap;

import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.io.template.FastaHeaderParserInterface;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

/**
 * Use FastaReaderHelper as an example of how to use this class where FastaReaderHelper should be the
 * primary class used to read Fasta files
 * @author Scooter Willis ;lt;willishf at gmail dot com&gt;
 */
public class FastaReader<S extends Sequence<?>, C extends Compound> {

    SequenceCreatorInterface<C> sequenceCreator;
    FastaHeaderParserInterface<S,C> headerParser;
    BufferedReaderBytesRead br;
    InputStreamReader isr;
    FileInputStream fi = null;
    long fileIndex = 0;
    long sequenceIndex = 0;
    String line = "";
    String header= "";
    
    /**
     * If you are going to use FileProxyProteinSequenceCreator then do not use this constructor because we need details about
     * local file offsets for quick reads. InputStreams does not give you the name of the stream to access quickly via file seek. A seek in
     * an inputstream is forced to read all the data so you don't gain anything.
     * @param br
     * @param headerParser
     * @param sequenceCreator
     */
    public FastaReader(InputStream is, FastaHeaderParserInterface<S,C> headerParser,
    		SequenceCreatorInterface<C> sequenceCreator) {
        this.headerParser = headerParser;
        isr = new InputStreamReader(is);
        this.br = new BufferedReaderBytesRead(isr);
        this.sequenceCreator = sequenceCreator;
    }

    /**
     * If you are going to use the FileProxyProteinSequenceCreator then you
     * need to use this constructor because we need details about
     * the location of the file.
     * @param file
     * @param headerParser
     * @param sequenceCreator
     * @throws FileNotFoundException if the file does not exist, is a directory 
     * 	rather than a regular file, or for some other reason cannot be opened
     * 	for reading.
     * @throws SecurityException if a security manager exists and its checkRead
     * 	method denies read access to the file.
     */
    public FastaReader(File file, FastaHeaderParserInterface<S,C> headerParser,
    		SequenceCreatorInterface<C> sequenceCreator) throws FileNotFoundException {
        this.headerParser = headerParser;
        fi = new FileInputStream(file);
        isr = new InputStreamReader(fi);
        this.br = new BufferedReaderBytesRead(isr);
        this.sequenceCreator = sequenceCreator;
    }

    /**
     * The parsing is done in this method.<br>
     * This method tries to process all the available fasta records 
     * in the File or InputStream, closes the underlying resource, 
     * and return the results in {@link LinkedHashMap}.<br>
     * You don't need to call {@link #close()} after calling this method.
     * @see #process(int)
     * @return {@link HashMap} containing all the parsed fasta records 
     * present, starting current fileIndex onwards.
     * @throws IOException if an error occurs reading the input file
     */
    public LinkedHashMap<String,S> process() throws IOException {
    	LinkedHashMap<String,S> sequences = process(-1);
    	close();
    	return sequences;
    }

    /**
     * This method tries to parse maximum <code>max</code> records from
     * the open File or InputStream, and leaves the underlying resource open.<br>
     * Subsequent calls to the same method continue parsing the rest of the file.<br>
     * This is particularly useful when dealing with very big data files,
     * (e.g. NCBI nr database), which can't fit into memory and will take long
     * time before the first result is available.<br>
     * <b>N.B.</b>
     * <ul>
     * <li>This method ca't be called after calling its NO-ARGUMENT twin.</li> 
     * <li>remember to close the underlying resource when you are done.</li> 
     * </ul>
     * @see #process()
     * @author Amr AL-Hossary
     * @since 3.0.6
     * @param max maximum number of records to return, <code>-1</code> for infinity.
     * @return {@link HashMap} containing maximum <code>max</code> parsed fasta records 
     * present, starting current fileIndex onwards.
     * @throws IOException if an error occurs reading the input file
     */
    public LinkedHashMap<String,S> process(int max) throws IOException {
        LinkedHashMap<String,S> sequences = new LinkedHashMap<String,S>();

        String line = "";
        if(this.line != null && this.line.length() > 0){
        	line=this.line;
        }
        String header = "";
        if(this.header != null && this.header.length() > 0){
        	header=this.header;
        }
        
        StringBuilder sb = new StringBuilder();
        int processedSequences=0;
        boolean keepGoing = true;

        do {
            line = line.trim(); // nice to have but probably not needed
            if (line.length() != 0) {
                if (line.startsWith(">")) {//start of new fasta record
                    if (sb.length() > 0) {//i.e. if there is already a sequence before
                    //    System.out.println("Sequence index=" + sequenceIndex);
                    	@SuppressWarnings("unchecked")
                        S sequence = (S)sequenceCreator.getSequence(sb.toString(), sequenceIndex);
                        headerParser.parseHeader(header, sequence);
                        sequences.put(sequence.getAccession().getID(),sequence);
                        processedSequences++;
//                        if (maxSequenceLength < sb.length()) {
//                            maxSequenceLength = sb.length();
//                        }
//                        sb = new StringBuilder(maxSequenceLength);
                        sb.setLength(0); //this is faster, better memory utilization (same buffer)
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
			if (line == null) {//i.e. EOF
				@SuppressWarnings("unchecked")
				//    System.out.println("Sequence index=" + sequenceIndex + " " + fileIndex );
                S sequence = (S)sequenceCreator.getSequence(sb.toString(), sequenceIndex);
                headerParser.parseHeader(header, sequence);
                sequences.put(sequence.getAccession().getID(),sequence);
                processedSequences++;
                keepGoing = false;
            }
			if (max > -1 && processedSequences>=max) {
				keepGoing=false;
			}
        } while (keepGoing);
        this.line  = line;
        this.header= header;
        return sequences;
    }

	public void close() throws IOException {
		br.close();
        isr.close();
        //If stream was created from File object then we need to close it
        if (fi != null) {
            fi.close();
        }
        this.line=this.header = null;
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
//                if(key.equals("Q98SJ1_CHICK/15-61")){
//                    int dummy = 1;
//                }
                System.out.println(proteinSequence.toString());

            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
