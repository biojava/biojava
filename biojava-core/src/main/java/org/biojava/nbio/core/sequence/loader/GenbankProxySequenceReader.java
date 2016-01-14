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
 * @author Karl Nicholas <github:karlnicholas>
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 08-08-2013
 *
 */
package org.biojava.nbio.core.sequence.loader;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.features.*;
import org.biojava.nbio.core.sequence.io.GenbankSequenceParser;
import org.biojava.nbio.core.sequence.io.GenericGenbankHeaderParser;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;

/**
 * @author Karl Nicholas <github:karlnicholas>
 * @author Jacek Grzebyta <github:jgrzebyta>
 */
public class GenbankProxySequenceReader<C extends Compound> extends StringProxySequenceReader<C> implements FeaturesKeyWordInterface, DatabaseReferenceInterface, FeatureRetriever {

	private final static Logger logger = LoggerFactory.getLogger(GenbankProxySequenceReader.class);

	private static final String eutilBaseURL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"; //
    private String genbankDirectoryCache = null;
    private GenbankSequenceParser<AbstractSequence<C>, C> genbankParser;
    private GenericGenbankHeaderParser<AbstractSequence<C>, C> headerParser;
    private String header;
    private HashMap<String, ArrayList<AbstractFeature>> features;
    

    /**
     * 
     * @throws InterruptedException 
     * @throws IOException 
     * @throws CompoundNotFoundException 
     */
    public GenbankProxySequenceReader(
            String genbankDirectoryCache,
            String accessionID,
            CompoundSet<C> compoundSet ) throws IOException, InterruptedException, CompoundNotFoundException {

        setGenbankDirectoryCache(genbankDirectoryCache);
        setCompoundSet(compoundSet);

        String db = compoundSet instanceof AminoAcidCompoundSet ? "protein" : "nuccore";

        InputStream inStream = getBufferedInputStream(accessionID, db);
        genbankParser = new GenbankSequenceParser<AbstractSequence<C>, C>();

        setContents(genbankParser.getSequence(new BufferedReader(new InputStreamReader(inStream)), 0));
        headerParser = genbankParser.getSequenceHeaderParser();
        header = genbankParser.getHeader();
        features = genbankParser.getFeatures();

        if (compoundSet.getClass().equals(AminoAcidCompoundSet.class)) {
            if (!genbankParser.getCompoundType().equals(compoundSet)) {
                logger.error("Declared compount type {} does not mach the real: {}", genbankParser.getCompoundType().toString(), compoundSet.toString());
                throw new IOException("Wrong declared compound type for: " + accessionID); 
            }
        }

        inStream.close();
    }

    private BufferedInputStream getBufferedInputStream(String accessionID, String db) throws IOException, InterruptedException {
        BufferedInputStream inStream = null;
        if (genbankDirectoryCache != null && genbankDirectoryCache.length() > 0) {
            File f = new File(genbankDirectoryCache + File.separatorChar + accessionID + ".gb");
            if (f.exists()) {
                logger.debug("Reading: {}", f.toString());
                inStream = new BufferedInputStream(new FileInputStream(f));
            } else {
                InputStream in = getEutilsInputStream(accessionID, db);
                copyInputStreamToFile(in, f);
                inStream = new BufferedInputStream(new FileInputStream(f));
            }
        } else {
            inStream = new BufferedInputStream(getEutilsInputStream(accessionID, db));
        }
        return inStream;
    }

    private void copyInputStreamToFile(InputStream in, File f) throws IOException, InterruptedException {
        FileOutputStream out = new FileOutputStream(f);
        byte[] buffer = new byte[1024];
        int len = in.read(buffer);
        while (len != -1) {
            out.write(buffer, 0, len);
            len = in.read(buffer);
            if (Thread.interrupted()) {
                in.close();
                out.close();
                throw new InterruptedException();
            }
        }
        in.close();
        out.close();
    }

    private InputStream getEutilsInputStream(String accessionID, String db) throws IOException {
        String genbankURL = eutilBaseURL + "efetch.fcgi?db=" + db + "&id=" + accessionID + "&rettype=gb&retmode=text";
        logger.trace("Loading: {}", genbankURL);
        URL genbank = new URL(genbankURL);
        URLConnection genbankConnection = genbank.openConnection();
        return genbankConnection.getInputStream();
    }

    /**
     * Local directory cache of Genbank that can be downloaded
     *
     * @return the uniprotDirectoryCache
     */
    public String getGenbankDirectoryCache() {
        return genbankDirectoryCache;
    }

    /**
     * @param genbankDirectoryCache
     */
    public void setGenbankDirectoryCache(String genbankDirectoryCache) {
        if (genbankDirectoryCache != null) {
            File f = new File(genbankDirectoryCache);
            if (!f.exists()) {
                f.mkdirs();
            }
        }
        this.genbankDirectoryCache = genbankDirectoryCache;
    }

    public String getHeader() {
        return header;
    }

    public GenericGenbankHeaderParser<AbstractSequence<C>, C> getHeaderParser() {
        return headerParser;
    }
    @Override
	public HashMap<String, ArrayList<AbstractFeature>> getFeatures() {
        return features;
    }

    @Override
    public LinkedHashMap<String, ArrayList<DBReferenceInfo>> getDatabaseReferences() {
        return genbankParser.getDatabaseReferences();
    }

    @Override
    public ArrayList<String> getKeyWords() {
        return genbankParser.getKeyWords();
    }

    public static void main(String[] args) throws Throwable {

        GenbankProxySequenceReader<AminoAcidCompound> genbankProteinReader
                = new GenbankProxySequenceReader<AminoAcidCompound>("/tmp", "NP_000257", AminoAcidCompoundSet.getAminoAcidCompoundSet());
        ProteinSequence proteinSequence = new ProteinSequence(genbankProteinReader);
        genbankProteinReader.getHeaderParser().parseHeader(genbankProteinReader.getHeader(), proteinSequence);
        logger.info("Sequence ({},{})={}...", proteinSequence.getAccession(), proteinSequence.getLength(), proteinSequence.getSequenceAsString().substring(0, 10));
        logger.info("Keywords: {}", genbankProteinReader.getKeyWords());
        logger.info("DatabaseReferences: {}", genbankProteinReader.getDatabaseReferences());
        proteinSequence.getFeatures();

        GenbankProxySequenceReader<NucleotideCompound> genbankDNAReader
                = new GenbankProxySequenceReader<NucleotideCompound>("/tmp", "NM_001126", DNACompoundSet.getDNACompoundSet());
        DNASequence dnaSequence = new DNASequence(genbankDNAReader);
        genbankDNAReader.getHeaderParser().parseHeader(genbankDNAReader.getHeader(), dnaSequence);
        dnaSequence.setAccession(new AccessionID("NM_001126"));
        logger.info("Sequence ({},{})={}...", dnaSequence.getAccession(), dnaSequence.getLength(), dnaSequence.getSequenceAsString().substring(0, 10));
        logger.info("Keywords: {}", genbankDNAReader.getKeyWords());
        logger.info("DatabaseReferences: {}", genbankDNAReader.getDatabaseReferences());

        genbankDNAReader
                = new GenbankProxySequenceReader<NucleotideCompound>("/tmp", "NM_000266", DNACompoundSet.getDNACompoundSet());
        dnaSequence = new DNASequence(genbankDNAReader);
        genbankDNAReader.getHeaderParser().parseHeader(genbankDNAReader.getHeader(), dnaSequence);
        logger.info("Sequence ({},{})={}...", dnaSequence.getAccession(), dnaSequence.getLength(), dnaSequence.getSequenceAsString().substring(0, 10));
        logger.info("Keywords: {}", genbankDNAReader.getKeyWords());
        logger.info("DatabaseReferences: {}", genbankDNAReader.getDatabaseReferences());

        genbankDNAReader
                = new GenbankProxySequenceReader<NucleotideCompound>("/tmp", "AV254721", DNACompoundSet.getDNACompoundSet());
        dnaSequence = new DNASequence(genbankDNAReader);
        genbankDNAReader.getHeaderParser().parseHeader(genbankDNAReader.getHeader(), dnaSequence);
        logger.info("Sequence ({},{})={}...", dnaSequence.getAccession(), dnaSequence.getLength(), dnaSequence.getSequenceAsString().substring(0, 10));
        logger.info("Keywords: {}", genbankDNAReader.getKeyWords());
        logger.info("DatabaseReferences: {}", genbankDNAReader.getDatabaseReferences());

        genbankDNAReader
                = new GenbankProxySequenceReader<NucleotideCompound>("/tmp", "AV254721.2", DNACompoundSet.getDNACompoundSet());
        dnaSequence = new DNASequence(genbankDNAReader);
        genbankDNAReader.getHeaderParser().parseHeader(genbankDNAReader.getHeader(), dnaSequence);
        logger.info("Sequence ({},{})={}...", dnaSequence.getAccession(), dnaSequence.getLength(), dnaSequence.getSequenceAsString().substring(0, 10));
        logger.info("Keywords: {}", genbankDNAReader.getKeyWords());
        logger.info("DatabaseReferences: {}", genbankDNAReader.getDatabaseReferences());

        genbankDNAReader
                = new GenbankProxySequenceReader<NucleotideCompound>("/tmp", "U49845", DNACompoundSet.getDNACompoundSet());
        dnaSequence = new DNASequence(genbankDNAReader);
        genbankDNAReader.getHeaderParser().parseHeader(genbankDNAReader.getHeader(), dnaSequence);
        logger.info("Sequence ({},{})={}...", dnaSequence.getAccession(), dnaSequence.getLength(), dnaSequence.getSequenceAsString().substring(0, 10));
        logger.info("Keywords: {}", genbankDNAReader.getKeyWords());
        logger.info("DatabaseReferences: {}", genbankDNAReader.getDatabaseReferences());

        genbankDNAReader
                = new GenbankProxySequenceReader<NucleotideCompound>("/tmp", "GI:1293613", DNACompoundSet.getDNACompoundSet());
        dnaSequence = new DNASequence(genbankDNAReader);
        genbankDNAReader.getHeaderParser().parseHeader(genbankDNAReader.getHeader(), dnaSequence);
        logger.info("Sequence ({},{})={}...", dnaSequence.getAccession(), dnaSequence.getLength(), dnaSequence.getSequenceAsString().substring(0, 10));
        logger.info("Keywords: {}", genbankDNAReader.getKeyWords());
        logger.info("DatabaseReferences: {}", genbankDNAReader.getDatabaseReferences());

        genbankDNAReader
                = new GenbankProxySequenceReader<NucleotideCompound>("/tmp", "14109166", DNACompoundSet.getDNACompoundSet());
        dnaSequence = new DNASequence(genbankDNAReader);
        genbankDNAReader.getHeaderParser().parseHeader(genbankDNAReader.getHeader(), dnaSequence);
        logger.info("Sequence ({},{})={}...", dnaSequence.getAccession(), dnaSequence.getLength(), dnaSequence.getSequenceAsString().substring(0, 10));
        logger.info("Keywords: {}", genbankDNAReader.getKeyWords());
        logger.info("DatabaseReferences: {}", genbankDNAReader.getDatabaseReferences());

        /*
         GenbankProxySequenceReader genbankProxyReader = new GenbankProxySequenceReader("/tmp");
         Sequence<?> sequence;

         sequence = genbankProxyReader.getDNASequence(new AccessionID("NM_001126"));
         System.out.println("Sequence" + "(" + sequence.getLength() + ")=" + sequence.getSequenceAsString().substring(0, 10) + "...");

         sequence = genbankProxyReader.getDNASequence(new AccessionID("NM_000266"));
         System.out.println("Sequence" + "(" + sequence.getLength() + ")=" + sequence.getSequenceAsString().substring(0, 10) + "...");
		
         sequence = genbankProxyReader.getProteinSequence(new AccessionID("NP_000257"));
         System.out.println("Sequence" + "(" + sequence.getLength() + ")=" + sequence.getSequenceAsString().substring(0, 10) + "...");
		
         sequence = genbankProxyReader.getProteinSequence(new AccessionID("AV254721"));
         System.out.println("Sequence" + "(" + sequence.getLength() + ")=" + sequence.getSequenceAsString().substring(0, 10) + "...");
		
         sequence = genbankProxyReader.getProteinSequence(new AccessionID("AV254721.2"));
         System.out.println("Sequence" + "(" + sequence.getLength() + ")=" + sequence.getSequenceAsString().substring(0, 10) + "...");
		
         sequence = genbankProxyReader.getProteinSequence(new AccessionID("U49845"));
         System.out.println("Sequence" + "(" + sequence.getLength() + ")=" + sequence.getSequenceAsString().substring(0, 10) + "...");
		
         sequence = genbankProxyReader.getProteinSequence(new AccessionID("GI:1293613"));
         System.out.println("Sequence" + "(" + sequence.getLength() + ")=" + sequence.getSequenceAsString().substring(0, 10) + "...");
		
         sequence = genbankProxyReader.getProteinSequence(new AccessionID("14109166"));
         System.out.println("Sequence" + "(" + sequence.getLength() + ")=" + sequence.getSequenceAsString().substring(0, 10) + "...");
         */
    }
}
