/**
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
 * Created on Dec 7, 2013
 * Created by Douglas Myers-Turnbull
 *
 */
package org.biojava.nbio.structure.io.sifts;

import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.core.sequence.io.util.IOUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.zip.GZIPInputStream;

/**
 * A mapping between UniProt entries and PDB chains.
 * For example
 * <pre>
 * SiftsChainToUniprot sifts = SiftsChainToUniprot.load();
 * SiftsChainEntry entry1 = sifts.getByUniProtId("P04585");
 * System.out.println(entry1.getPdbId() + "." + entry1.getChainName()); // 1hiv.A
 * System.out.println(entry1.getPdbStart() + "-" + entry1.getPdbStop()); // 1-99
 * SiftsChainEntry entry2 = sifts.getByChainId("1hiv", "A");
 * System.out.println(entry1.equals(entry2)); // true
 * </pre>
 * See SIFTS project documentation: https://www.ebi.ac.uk/pdbe/docs/sifts/
 * @author dmyersturnbull
 * @see SiftsChainEntry
 * @since 3.0.7
 */
public class SiftsChainToUniprotMapping {

	private final static Logger logger = LoggerFactory.getLogger(SiftsChainToUniprotMapping.class);


	protected static File DEFAULT_FILE;

	private static final String DEFAULT_FILENAME = "pdb_chain_uniprot.tsv";
	private static final URL DEFAULT_URL;

	static {
		try {
			DEFAULT_URL = new URL("http://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz");
		} catch (MalformedURLException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Loads the SIFTS mapping.
	 * Attempts to load the mapping file in the PDB cache directory.
	 * If the file does not exist or could not be parsed, downloads and stores a GZ-compressed file.
	 * @return
	 * @throws IOException If the local file could not be read and could not be downloaded
	 */
	public static SiftsChainToUniprotMapping load() throws IOException {
		return load(false);
	}

	/**
	 * Loads the SIFTS mapping.
	 * Attempts to load the mapping file in the PDB cache directory.
	 * If the file does not exist or could not be parsed, downloads and stores a GZ-compressed file.
	 * @param useOnlyLocal If true, will throw an IOException if the file needs to be downloaded
	 * @return
	 * @throws IOException If the local file could not be read and could not be downloaded (including if onlyLocal is true)
	 */
	public static SiftsChainToUniprotMapping load(boolean useOnlyLocal) throws IOException {

		UserConfiguration config = new UserConfiguration();
		File cacheDir = new File(config.getCacheFilePath());

		DEFAULT_FILE = new File(cacheDir, DEFAULT_FILENAME);


		if (!DEFAULT_FILE.exists() || DEFAULT_FILE.length() == 0) {
			if (useOnlyLocal) throw new IOException(DEFAULT_FILE + " does not exist, and did not download");
			download();
		}
		try {
			return build();
		} catch (IOException e) {
			logger.info("Caught IOException while reading {}. Error: {}",DEFAULT_FILE,e.getMessage());
			if (useOnlyLocal) throw new IOException(DEFAULT_FILE + " could not be read, and did not redownload");
			download();
			return build();
		}
	}

	/**
	 * Builds the mapping by reading SIFTS the tsv file set in {@link #DEFAULT_FILE} variable.
	 * @return
	 * @throws IOException
	 */
	protected static SiftsChainToUniprotMapping build() throws IOException {
		SiftsChainToUniprotMapping sifts = new SiftsChainToUniprotMapping();
		BufferedReader br = new BufferedReader(new FileReader(DEFAULT_FILE));
		String line = "";
		while ((line = br.readLine()) != null) {
			if (line.isEmpty() || line.startsWith("#") || line.startsWith("PDB")) continue;
			String[] parts = line.split("\t");
			String pdbId = parts[0];
			String chainId = parts[1];
			String uniProtId = parts[2];
			String seqresStart = parts[3];
			String seqresEnd = parts[4];
			String pdbStart = parts[5];
			String pdbEnd = parts[6];
			String uniprotStart = parts[7];
			String uniprotEnd = parts[8];
			SiftsChainEntry entry = new SiftsChainEntry(pdbId, chainId, uniProtId, seqresStart, seqresEnd,
					pdbStart, pdbEnd, uniprotStart, uniprotEnd);
			sifts.byChainId.put(pdbId + "." + chainId, entry);
			sifts.byUniProtId.put(uniProtId, entry);
		}
		br.close();
		return sifts;
	}

	private static void download() throws IOException {

		logger.info("Downloading {} to {}",DEFAULT_URL.toString(),DEFAULT_FILE);

		InputStream in = null;
		OutputStream out = null;

		in = new GZIPInputStream(DEFAULT_URL.openStream());
		out = new FileOutputStream(DEFAULT_FILE);
		IOUtils.copy(in, out);

	}

	private Map<String, SiftsChainEntry> byChainId = new HashMap<String, SiftsChainEntry>();

	private Map<String, SiftsChainEntry> byUniProtId = new HashMap<String, SiftsChainEntry>();

	private SiftsChainToUniprotMapping() {

	}

	public Set<Entry<String, SiftsChainEntry>> chainEntrySet() {
		return byChainId.entrySet();
	}

	public boolean containsChainId(String pdbId, String chainId) {
		return byChainId.containsKey(pdbId + "." + chainId);
	}

	public boolean containsUniProtId(String uniProtId) {
		return byUniProtId.containsKey(uniProtId);
	}

	public SiftsChainEntry getByChainId(String pdbId, String chainId) {
		return byChainId.get(pdbId + "." + chainId);
	}

	public SiftsChainEntry getByUniProtId(String uniProtId) {
		return byUniProtId.get(uniProtId);
	}

	public Set<String> keySet() {
		return byChainId.keySet();
	}

	/**
	 * Returns the number of mapped entries.
	 */
	public int size() {
		return byChainId.size();
	}

	public Set<Entry<String, SiftsChainEntry>> uniProtEntrySet() {
		return byChainId.entrySet();
	}

	public Collection<SiftsChainEntry> values() {
		return byChainId.values();
	}
}
