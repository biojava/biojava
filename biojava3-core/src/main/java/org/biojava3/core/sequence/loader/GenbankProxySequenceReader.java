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
 * Created on 08-08-2013
 *
 * @author Karl Nicholas
 *
 */
package org.biojava3.core.sequence.loader;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.net.URLConnection;
import java.util.LinkedHashMap;
import java.util.logging.Logger;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.ProteinSequence;

import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.DNASequenceCreator;
import org.biojava3.core.sequence.io.GenbankReader;
import org.biojava3.core.sequence.io.GenericGenbankHeaderParser;
import org.biojava3.core.sequence.io.ProteinSequenceCreator;

/**
 * 
 */
public class GenbankProxySequenceReader {

	private static final Logger logger = Logger.getLogger(UniprotProxySequenceReader.class.getName());
	private static final String eutilBaseURL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"; //
	private String genbankDirectoryCache = null;

	/**
	 * 
	 */
	public GenbankProxySequenceReader() {
	}

	/**
	 * 
	 */
	public GenbankProxySequenceReader(String genbankDirectoryCache) {
		setGenbankDirectoryCache(genbankDirectoryCache);
	}

	/**
	 * 
	 * Pass in a AccessionID and return a 
	 * DNASequence will get the sequence data and other data elements associated
	 * with the DNASequence in the Genbank format. This is an example of how to map
	 * external databases of nucleotides and features to the BioJava3 DNASequence.
	 * Important to call @see setGenebankDirectoryCache to allow caching of XML files
	 * so they don't need to be reloaded each time. Does not manage cache.
	 * 
	 * @param AccessionID
	 * @return ProteinSequence
	 * @throws IOException, InterruptedException
	 */
	public DNASequence getDNASequence(AccessionID accessionID) throws IOException, InterruptedException {
		BufferedInputStream inStream = getBufferedInputStream(accessionID, "nuccore" );
		GenbankReader<DNASequence, NucleotideCompound> genbankReader = new GenbankReader<DNASequence, NucleotideCompound>(
				inStream,
				new GenericGenbankHeaderParser<DNASequence, NucleotideCompound>(),
				new DNASequenceCreator(DNACompoundSet.getDNACompoundSet()));
		LinkedHashMap<String, DNASequence> dnaSequences = genbankReader.process();
		inStream.close();
		return dnaSequences.get(accessionID.getID());
	}

	/**
	 * 
	 * Pass in a AccessionID and return a 
	 * ProteinSequence will get the sequence data and other data elements associated
	 * with the ProteinSequence in the Genbank format. This is an example of how to map
	 * external databases of proteins and features to the BioJava3 ProteinSequence.
	 * Important to call @see setGenebankDirectoryCache to allow caching of Genbank files
	 * so they don't need to be reloaded each time. Does not manage cache.
	 * 
	 * @param AccessionID
	 * @return ProteinSequence
	 * @throws IOException, InterruptedException
	 */
	public ProteinSequence getProteinSequence(AccessionID accessionID) throws IOException, InterruptedException {
		BufferedInputStream inStream = getBufferedInputStream(accessionID, "protein" );
		GenbankReader<ProteinSequence, AminoAcidCompound> GenbankProtein = new GenbankReader<ProteinSequence, AminoAcidCompound>(
				inStream,
				new GenericGenbankHeaderParser<ProteinSequence, AminoAcidCompound>(),
				new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
		LinkedHashMap<String, ProteinSequence> proteinSequences = GenbankProtein.process();
		inStream.close();
		return proteinSequences.get(accessionID.getID());
	}
	
	private BufferedInputStream getBufferedInputStream(AccessionID accessionID, String db) throws IOException, InterruptedException {
		BufferedInputStream inStream = null;
		if (genbankDirectoryCache != null && genbankDirectoryCache.length() > 0) {
			File f = new File(genbankDirectoryCache + File.separatorChar + accessionID + ".gb");
			if (f.exists()) {
				logger.info("Reading " + f.toString());
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

	private InputStream getEutilsInputStream(AccessionID accessionID, String db) throws IOException {
		String genbankURL = eutilBaseURL + "efetch.fcgi?db=" + db + "&id=" + accessionID + "&rettype=gb&retmode=text";
		logger.info("Loading " + genbankURL);
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
	 * @param aGenbankDirectoryCache the GenbankDirectoryCache to set
	 */
	public void setGenbankDirectoryCache(String genbankDirectoryCache) {
		File f = new File(genbankDirectoryCache);
		if (!f.exists()) {
			f.mkdirs();
		}
		this.genbankDirectoryCache = genbankDirectoryCache;
	}

	public static void main(String[] args) throws Exception {
		GenbankProxySequenceReader genbankProxyReader = new GenbankProxySequenceReader("/tmp/dna");

		DNASequence dnaSequence = genbankProxyReader.getDNASequence(new AccessionID("NM_001126"));
		System.out.println("Sequence=" + dnaSequence.getSequenceAsString());
		dnaSequence = genbankProxyReader.getDNASequence(new AccessionID("NM_000266"));
		System.out.println("Sequence=" + dnaSequence.getSequenceAsString());

		genbankProxyReader.setGenbankDirectoryCache("/tmp/aa");
		ProteinSequence proteinSequence = genbankProxyReader.getProteinSequence(new AccessionID("NP_000257"));
		System.out.println("Sequence=" + proteinSequence.getSequenceAsString());
	}

}