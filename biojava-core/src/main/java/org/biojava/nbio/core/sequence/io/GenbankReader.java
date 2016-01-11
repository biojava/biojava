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
 * @author Scooter Willis ;lt;willishf at gmail dot com&gt;
 * @author Karl Nicholas <github:karlnicholas>
 * @author Paolo Pavan
 * 
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 01-21-2010
 */
package org.biojava.nbio.core.sequence.io;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.DataSource;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.TaxonomyID;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.features.AbstractFeature;
import org.biojava.nbio.core.sequence.features.DBReferenceInfo;
import org.biojava.nbio.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava.nbio.core.sequence.io.template.SequenceHeaderParserInterface;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.template.Compound;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;

/**
 * Use GenbankReaderHelper as an example of how to use this class where GenbankReaderHelper should be the
 * primary class used to read Genbank files
 *
 */
public class GenbankReader<S extends AbstractSequence<C>, C extends Compound> {

    private SequenceCreatorInterface<C> sequenceCreator;
    private GenbankSequenceParser<S,C> genbankParser;
    private InputStream inputStream;
    
    /**
     * If you are going to use FileProxyProteinSequenceCreator then do not use this constructor because we need details about
     * local file offsets for quick reads. InputStreams does not give you the name of the stream to access quickly via file seek. A seek in
     * an inputstream is forced to read all the data so you don't gain anything.
     * @param br
     * @param headerParser
     * @param sequenceCreator
     */
    public GenbankReader(InputStream is, SequenceHeaderParserInterface<S,C> headerParser, SequenceCreatorInterface<C> sequenceCreator) {
        this.sequenceCreator = sequenceCreator;
        this.inputStream = is;
    	genbankParser = new GenbankSequenceParser<S,C>();
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
    public GenbankReader(
    		File file, 
    		SequenceHeaderParserInterface<S,C> headerParser, 
    		SequenceCreatorInterface<C> sequenceCreator
    		) throws FileNotFoundException {
    	
        inputStream = new FileInputStream(file);
        this.sequenceCreator = sequenceCreator;
    	genbankParser = new GenbankSequenceParser<S,C>();
    }

    /**
     * The parsing is done in this method.<br>
     * This method tries to process all the available Genbank records 
     * in the File or InputStream, closes the underlying resource, 
     * and return the results in {@link LinkedHashMap}.<br>
     * You don't need to call {@link #close()} after calling this method.
     * @see #process(int)
     * @return {@link HashMap} containing all the parsed Genbank records 
     * present, starting current fileIndex onwards.
     * @throws IOException
     * @throws CompoundNotFoundException 
     */
    public LinkedHashMap<String,S> process() throws IOException, CompoundNotFoundException {
    	LinkedHashMap<String,S> sequences = process(-1);
    	return sequences;
    }

    /**
     * This method tries to parse maximum <code>max</code> records from
     * the open File or InputStream, and leaves the underlying resource open.<br>
     * 
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
     * @return {@link HashMap} containing maximum <code>max</code> parsed Genbank records 
     * present, starting current fileIndex onwards.
     * @throws IOException
     * @throws CompoundNotFoundException 
     */
    public LinkedHashMap<String,S> process(int max) throws IOException, CompoundNotFoundException {
        LinkedHashMap<String,S> sequences = new LinkedHashMap<String,S>();
        @SuppressWarnings("unchecked")
        int i=0;
        BufferedReader br = new BufferedReader(new InputStreamReader(inputStream));
        while(true) {
        	if(max>0 && i>=max) break;
        	i++;
        	String seqString = genbankParser.getSequence(br, 0);
        	//reached end of file?
        	if(seqString==null) break;
        	S sequence = (S) sequenceCreator.getSequence(seqString, 0);
        	genbankParser.getSequenceHeaderParser().parseHeader(genbankParser.getHeader(), sequence);
        	
        	// add features to new sequence
        	for (String k: genbankParser.getFeatures().keySet()){
        		for (AbstractFeature f: genbankParser.getFeatures(k)){
        			//f.getLocations().setSequence(sequence);  // can't set proper sequence source to features. It is actually needed? Don't think so...
        			sequence.addFeature(f);
        		}
        	}
        
        	// add taxonomy ID to new sequence
        	ArrayList<DBReferenceInfo> dbQualifier = genbankParser.getDatabaseReferences().get("db_xref");
        	if (dbQualifier != null){
        		DBReferenceInfo q = dbQualifier.get(0);
        		sequence.setTaxonomy(new TaxonomyID(q.getDatabase()+":"+q.getId(), DataSource.GENBANK));
        	}
        	
        	sequences.put(sequence.getAccession().getID(), sequence);
        }
    	br.close();
    	close();
        return sequences;
    }

	public void close() throws IOException {
		inputStream.close();
	}

    public static void main(String[] args) throws Exception {
        String proteinFile = "src/test/resources/BondFeature.gb";
        FileInputStream is = new FileInputStream(proteinFile);

        GenbankReader<ProteinSequence, AminoAcidCompound> proteinReader = new GenbankReader<ProteinSequence, AminoAcidCompound>(is, new GenericGenbankHeaderParser<ProteinSequence,AminoAcidCompound>(), new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
        LinkedHashMap<String,ProteinSequence> proteinSequences = proteinReader.process();
        System.out.println(proteinSequences);

        String inputFile = "src/test/resources/NM_000266.gb";
        is = new FileInputStream(inputFile);
        GenbankReader<DNASequence, NucleotideCompound> dnaReader = new GenbankReader<DNASequence, NucleotideCompound>(is, new GenericGenbankHeaderParser<DNASequence,NucleotideCompound>(), new DNASequenceCreator(DNACompoundSet.getDNACompoundSet()));
        LinkedHashMap<String,DNASequence> dnaSequences = dnaReader.process();
        System.out.println(dnaSequences);
        
        String crazyFile = "src/test/resources/CraftedFeature.gb";
        is = new FileInputStream(crazyFile);
        GenbankReader<DNASequence, NucleotideCompound> crazyReader = new GenbankReader<DNASequence, NucleotideCompound>(is, new GenericGenbankHeaderParser<DNASequence,NucleotideCompound>(), new DNASequenceCreator(DNACompoundSet.getDNACompoundSet()));
        LinkedHashMap<String,DNASequence> crazyAnnotatedSequences = crazyReader.process();
        
        is.close();
        System.out.println(crazyAnnotatedSequences);
    }

}

