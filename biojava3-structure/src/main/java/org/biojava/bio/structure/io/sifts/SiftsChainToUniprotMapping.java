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
package org.biojava.bio.structure.io.sifts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava3.core.sequence.io.util.IOUtils;

/**
 * A mapping between UniProt entries and PDB chains.
 * For example
 * <pre>
 * SiftsChainToUniprot sifts = SiftsChainToUniprot.load();
 * SiftsChainEntry entry1 = sifts.getByUniProtId("P04585");
 * System.out.println(entry1.getPdbId() + "." + entry1.getChainId()); // 1hiv.A
 * System.out.println(entry1.getPdbStart() + "-" + entry1.getPdbStop()); // 1-99
 * SiftsChainEntry entry2 = sifts.getByChainId("1hiv", "A");
 * System.out.println(entry1.equals(entry2)); // true
 * </pre>
 * 
 * @author dmyersturnbull
 * @see SiftsChainEntry
 * @since 3.0.7
 */
public class SiftsChainToUniprotMapping {

	private static File DEFAULT_FILE;

	private static final String DEFAULT_FILENAME = "pdb_chain_uniprot.tsv";
	private static final URL DEFAULT_URL;

	static {
		try {
			DEFAULT_URL = new URL("ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz");
		} catch (MalformedURLException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Loads the SIFTS mapping.
	 * Attempts to load the mapping file file in the PDB cache directory.
	 * If the file does not exist or could not be parsed, downloads and stores a GZ-compressed file.
	 * @return
	 * @throws IOException If the local file could not be read and could not be downloaded
	 */
	public static SiftsChainToUniprotMapping load() throws IOException {
		return load(false);
	}

	/**
	 * Loads the SIFTS mapping.
	 * Attempts to load the mapping file file in the PDB cache directory.
	 * If the file does not exist or could not be parsed, downloads and stores a GZ-compressed file.
	 * @param useOnlyLocal If true, will throw an IOException if the file needs to be downloaded
	 * @return
	 * @throws IOException If the local file could not be read and could not be downloaded (including if onlyLocal is true)
	 */
	public static SiftsChainToUniprotMapping load(boolean useOnlyLocal) throws IOException {
		String cacheDir = System.getProperty(AbstractUserArgumentProcessor.CACHE_DIR);
		if (cacheDir != null) {
			DEFAULT_FILE = new File(cacheDir.endsWith(File.pathSeparator) ? cacheDir + DEFAULT_FILENAME : cacheDir
					+ File.pathSeparator + DEFAULT_FILENAME);
		} else {
			DEFAULT_FILE = File.createTempFile(DEFAULT_FILENAME, "xml");
		}
		if (!DEFAULT_FILE.exists()) {
			if (useOnlyLocal) throw new IOException(DEFAULT_FILE + " does not exist, and did not download");
			download();
		}
		try {
			return build();
		} catch (IOException e) {
			e.printStackTrace();
			if (useOnlyLocal) throw new IOException(DEFAULT_FILE + " could not be read, and did not redownload");
			download();
			return build();
		}
	}

	private static SiftsChainToUniprotMapping build() throws IOException {
		SiftsChainToUniprotMapping sifts = new SiftsChainToUniprotMapping();
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(DEFAULT_FILE));
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
		} finally {
			if (br != null) {
				try {
					br.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		return sifts;
	}

	private static void download() throws IOException {
		InputStream in = null;
		OutputStream out = null;
		try {
			in = new GZIPInputStream(DEFAULT_URL.openStream());
			out = new FileOutputStream(DEFAULT_FILE);
			IOUtils.copy(in, out);
		} finally {
			if (in != null) {
				try {
					in.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			if (in != null) {
				try {
					out.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
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
