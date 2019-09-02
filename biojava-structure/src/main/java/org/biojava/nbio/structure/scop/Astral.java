/**
 * BioJava development code
 *
 * This code may be freely distributed and modified under the terms of the GNU Lesser General Public Licence. This
 * should be distributed with the code. If you do not have a copy, see:
 *
 * http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual authors. These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims, or to join the biojava-l mailing list, visit the home page
 * at:
 *
 * http://www.biojava.org/
 *
 */
package org.biojava.nbio.structure.scop;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.lang.ref.SoftReference;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;


/**
 * Provides programmatic access to ASTRAL representative sets. See the paper by <a
 * href="http://scop.berkeley.edu/references/2004-nar-astral.pdf">Chandonia et. al.</a> for more information. Example:
 *
 * <pre>
 * Set&lt;String&gt; astralSet = Astral.getRepresentatives(Astral.AstralSet.NINETY_FIVE_175B);
 * </pre>
 *
 * This class uses a multiton pattern with soft references for caching. In short: the first time you call the above, it
 * will fetch the data from ASTRAL; the second time will (probably) not have to; and the instances can still be
 * garbage-collected if necessary (meaning they don't <em>require</em> heap memory).
 *
 * @author dmyerstu
 * @since 3.0.6
 */
public class Astral {

	/**
	 * An ASTRAL sequence-identity cutoff with an identifier such as:
	 *
	 * <pre>
	 * 1.75A_95
	 * </pre>
	 *
	 * Also contains a URL pointing to a FASTA file containing the representatives. Every character before the first
	 * whitespace character of each header in the FASTA file is expected to be a representative's name.
	 *
	 * @author dmyersturnbull
	 *
	 */
	public static enum AstralSet {
		FORTY_175("1.75_40", "http://scop.berkeley.edu/downloads/scopseq-1.75/astral-scopdom-seqres-gd-sel-gs-bib-40-1.75.fa"),
		NINETY_FIVE_175("1.75_95", "http://scop.berkeley.edu/downloads/scopseq-1.75/astral-scopdom-seqres-gd-sel-gs-bib-95-1.75.fa"),
		FORTY_175A("1.75A_40", "http://scop.berkeley.edu/downloads/scopeseq-2.01/astral-scopedom-seqres-gd-sel-gs-bib-40-2.01.fa"),
		NINETY_FIVE_175A("1.75A_95","http://scop.berkeley.edu/downloads/scopeseq-2.01/astral-scopedom-seqres-gd-sel-gs-bib-95-2.01.fa"),
		FORTY_175B("1.75B_40", "http://scop.berkeley.edu/downloads/scopeseq-2.02/astral-scopedom-seqres-gd-sel-gs-bib-40-2.02.fa"),
		NINETY_FIVE_175B("1.75B_95", "http://scop.berkeley.edu/downloads/scopeseq-2.02/astral-scopedom-seqres-gd-sel-gs-bib-95-2.02.fa"),
		FORTY_201("2.01_40", "http://scop.berkeley.edu/downloads/scopeseq-2.01/astral-scopedom-seqres-gd-sel-gs-bib-40-2.01.fa"),
		NINETY_FIVE_201("2.01_95", "http://scop.berkeley.edu/downloads/scopeseq-2.01/astral-scopedom-seqres-gd-sel-gs-bib-95-2.01.fa"),
		FORTY_202("2.02_40", "http://scop.berkeley.edu/downloads/scopeseq-2.02/astral-scopedom-seqres-gd-sel-gs-bib-40-2.02.fa"),
		NINETY_FIVE_202("2.02_95", "http://scop.berkeley.edu/downloads/scopeseq-2.02/astral-scopedom-seqres-gd-sel-gs-bib-95-2.02.fa"),
		FORTY_203("2.03_40", "http://scop.berkeley.edu/downloads/scopeseq-2.03/astral-scopedom-seqres-gd-sel-gs-bib-40-2.03.fa"),
		NINETY_FIVE_203("2.03_95", "http://scop.berkeley.edu/downloads/scopeseq-2.03/astral-scopedom-seqres-gd-sel-gs-bib-95-2.03.fa");
		private String id;
		private String url;

		public static AstralSet parse(String str) {
			for (AstralSet c : AstralSet.class.getEnumConstants()) {
				if (c.getId().equals(str)) return c;
			}
			throw new IllegalArgumentException("No ASTRAL set with id " + str);
		}

		AstralSet(String id, String url) {
			this.url = url;
			this.id = id;
		}

		public String getId() {
			return id;
		}

		public String getUrl() {
			return url;
		}

		@Override
		public String toString() {
			return id;
		}
	}

	private static Map<String, SoftReference<Astral>> instances = new HashMap<String, SoftReference<Astral>>();

	private static final Logger logger = LoggerFactory.getLogger(Astral.class);

	private Set<String> names;
	private LinkedHashMap<Integer,String> failedLines;

	/**
	 * Get a list of representatives' names for the specified ASTRAL cutoff.
	 */
	public static Set<String> getRepresentatives(AstralSet cutoff) {
		if (instances.containsKey(cutoff.getId()) && instances.get(cutoff.getId()).get() != null) {
			return instances.get(cutoff.getId()).get().getNames();
		}
		Astral astral = new Astral(cutoff);
		instances.put(cutoff.getId(), new SoftReference<Astral>(astral));
		return astral.getNames();
	}

	/**
	 * Get a list of representatives' names for the specified ASTRAL cutoff.
	 * @param id An ASTRAL Id, such as 1.75A_95.
	 */
	public static Set<String> getRepresentatives(String id) {
		return getRepresentatives(AstralSet.parse(id));
	}

	/**
	 * Constructs a new Astral object. Generally, client code should prefer calling
	 * {@link #getRepresentatives(AstralSet)} instead. This constructor should only be used when an ASTRAL set not
	 * included in {@link #Astral(AstralSet)} is required.
	 *
	 * @param cutoff
	 *            The ASTRAL sequence-identity cutoff required
	 * @throws RuntimeException
	 *             If the Astral set could not be parsed or accessed for any reason
	 */
	public Astral(AstralSet cutoff) {
		URL url;
		try {
			url = new URL(cutoff.getUrl());
		} catch (MalformedURLException e) {
			throw new RuntimeException("The URL was invalid!", e);
		}
		Reader reader;
		try {
			reader = new InputStreamReader(url.openStream());
		} catch (IOException e) {
			throw new RuntimeException("Couldn't open stream to URL " + url, e);
		}
		init(reader);
	}

	/**
	 * Constructs a new Astral object. Generally, client code should prefer calling
	 * {@link #getRepresentatives(AstralSet)} instead. This constructor should only be used when an ASTRAL set not
	 * included in {@link #Astral(AstralSet)} is required.
	 *
	 * @throws RuntimeException
	 *             If the Astral set could not be parsed or accessed for any reason
	 */
	public Astral(String id, URL url) {
		Reader reader;
		try {
			reader = new InputStreamReader(url.openStream());
		} catch (IOException e) {
			throw new RuntimeException("Couldn't open stream to URL " + url, e);
		}
		init(reader);
	}

	/**
	 * Constructs a new Astral object. Generally, client code should prefer calling
	 * {@link #getRepresentatives(AstralSet)} instead. This constructor should only be used when an ASTRAL set not
	 * included in {@link #Astral(AstralSet)} is required.
	 *
	 * @throws RuntimeException
	 *             If the Astral set could not be parsed or accessed for any reason
	 */
	public Astral(String id, Reader reader) {
		init(reader);
	}

	/**
	 * @return The names of representatives in this ASTRAL set.
	 */
	public Set<String> getNames() {
		return names;
	}

	/**
	 * Gets a map describing lines read in the file that weren't understood.
	 * @return A LinkedHashMap mapping line numbers of failures to the lines themselves
	 */
	public LinkedHashMap<Integer, String> getFailedLines() {
		return failedLines;
	}

	/**
	 * Parses the FASTA file opened by reader.
	 */
	private void init(Reader reader) {
		names = new TreeSet<String>();
		failedLines = new LinkedHashMap<Integer,String>();

		BufferedReader br = null;

		try {

			br = new BufferedReader(reader);

			logger.info("Reading ASTRAL file...");

			String line = "";
			int i = 0;
			while ((line = br.readLine()) != null) {
				if (line.startsWith(">")) {
					try {
						String scopId = line.split("\\s")[0].substring(1);
						names.add(scopId);
						if (i % 1000 == 0) {
							logger.debug("Reading ASTRAL line for " + scopId);
						}
						i++;
					} catch (RuntimeException e) {
						failedLines.put(i, line);
						logger.warn("Couldn't read line " + line, e);
					}
				}
			}

			br.close();

		} catch (IOException e) {
			throw new RuntimeException("Couldn't read the input stream ", e);
		} finally {
			if (br != null) {
				try {
					br.close();
				} catch (IOException e) {
					logger.warn("Could not close stream", e);
				}
			}
		}

	}

}
