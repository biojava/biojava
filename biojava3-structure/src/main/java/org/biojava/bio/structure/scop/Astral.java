/**
 *                  BioJava development code
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
 */
package org.biojava.bio.structure.scop;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Set;
import java.util.TreeSet;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Provides programmatic access to ASTRAL representative sets.
 * See the paper by <a href="http://scop.berkeley.edu/references/2004-nar-astral.pdf">Chandonia et. al.</a> for more information.
 * Example:
 * <pre>
 * Set&lt;String&gt; astralSet = Astral.getRepresentatives(Astral.AstralSet.NINETY_FIVE_175B);
 * </pre>
 * @author dmyerstu
 * @since 3.0.6
 */
public class Astral {

	private static final Logger logger = LogManager.getLogger(Astral.class.getPackage().getName());

	public static enum AstralSet {
		FORTY_175A("1.75A_40", "http://scop.berkeley.edu/downloads/scopseq-1.75A/astral-scopdom-seqres-gd-sel-gs-bib-40-1.75A.fa"),
		NINETY_FIVE_175A("1.75A_95", "http://scop.berkeley.edu/downloads/scopseq-1.75A/astral-scopdom-seqres-gd-sel-gs-bib-90-1.75A.fa"),
		FORTY_175B("1.75B_95", "http://scop.berkeley.edu/downloads/scopseq-1.75B/astral-scopdom-seqres-gd-sel-gs-bib-40-1.75B.fa"),
		NINETY_FIVE_175B("1.75B_95", "http://scop.berkeley.edu/downloads/scopseq-1.75B/astral-scopdom-seqres-gd-sel-gs-bib-95-1.75B.fa"),
		FORTY_175("1.75_95", "http://scop.berkeley.edu/downloads/scopseq-1.75/astral-scopdom-seqres-gd-sel-gs-bib-40-1.75.fa"),
		NINETY_FIVE_175("1.75_95", "http://scop.berkeley.edu/downloads/scopseq-1.75/astral-scopdom-seqres-gd-sel-gs-bib-95-1.75.fa");
		private String url;
		private String id;
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
		String getId() {
			return id;
		}
		String getUrl() {
			return url;
		}
	}

	public static Set<String> getRepresentatives(AstralSet cutoff) {
		
		URL url;
		try {
			url = new URL(cutoff.getUrl());
		} catch (MalformedURLException e) {
			throw new RuntimeException("The URL was invalid!", e);
		}
		
		Set<String> names = new TreeSet<String>();
		
		try {
			
			BufferedReader br = new BufferedReader(new InputStreamReader(url.openStream()));

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
						logger.error("Couldn't read line " + line, e);
					}
				}
			}

			br.close();
			
		} catch (IOException e) {
			throw new RuntimeException("Couldn't read the input stream to " + url.getPath(), e);
		}

		return names;
	}

}
