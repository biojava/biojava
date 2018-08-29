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
 */
package org.biojava.nbio.structure.symmetry.utils;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Serializable;
import java.net.URL;
import java.util.*;


public class BlastClustReader implements Serializable {

	private static final long serialVersionUID = 1L;

	private static final Logger logger = LoggerFactory.getLogger(BlastClustReader.class);

	private int sequenceIdentity = 0;
	private List<List<String>> clusters = new ArrayList<>();
	// https://cdn.rcsb.org/resources/sequence/clusters/bc-95.out
	private static final String coreUrl = "https://cdn.rcsb.org/resources/sequence/clusters/";

	private static final List<Integer> seqIdentities = Arrays.asList(30, 40, 50, 70, 90, 95, 100);

	public BlastClustReader(int sequenceIdentity)  {
		this.sequenceIdentity = sequenceIdentity;
	}

	public List<List<String>> getPdbChainIdClusters() {
		loadClusters(sequenceIdentity);
		return clusters;
	}

	public Map<String,String> getRepresentatives(String pdbId) {
		loadClusters(sequenceIdentity);
		String pdbIdUc = pdbId.toUpperCase();

		Map<String,String> representatives = new LinkedHashMap<>();
		for (List<String> cluster: clusters) {
			// map fist match to representative
			for (String chainId: cluster) {
				if (chainId.startsWith(pdbIdUc)) {
					representatives.put(chainId, cluster.get(0));
					break;
				}
			}
		}
		return representatives;
	}

	public String getRepresentativeChain(String pdbId, String chainId) {
		loadClusters(sequenceIdentity);

		String pdbChainId = pdbId.toUpperCase() + "." + chainId;

		for (List<String> cluster: clusters) {
			if (cluster.contains(pdbChainId)) {
				return cluster.get(0);
			}
		}
		return "";
	}

	public int indexOf(String pdbId, String chainId) {
		loadClusters(sequenceIdentity);

		String pdbChainId = pdbId.toUpperCase() + "." + chainId;

		for (int i = 0; i < clusters.size(); i++) {
			List<String> cluster = clusters.get(i);
			if (cluster.contains(pdbChainId)) {
				return i;
			}
		}
		return -1;
	}

	public List<List<String>> getPdbChainIdClusters(String pdbId) {
		loadClusters(sequenceIdentity);
		String pdbIdUpper = pdbId.toUpperCase();

		List<List<String>> matches = new ArrayList<List<String>>();
		for (List<String> cluster: clusters) {
			for (String chainId: cluster) {
				if (chainId.startsWith(pdbIdUpper)) {
					matches.add(cluster);
					break;
				}
			}
		}
		return matches;
	}

	public List<List<String>> getChainIdsInEntry(String pdbId) {
		loadClusters(sequenceIdentity);

		List<List<String>> matches = new ArrayList<List<String>>();
		List<String> match = null;

		for (List<String> cluster: clusters) {
			for (String chainId: cluster) {
				if (chainId.startsWith(pdbId)) {
					if (match == null) {
						match = new ArrayList<String>();
					}
					match.add(chainId.substring(5));
				}
			}
			if (match != null) {
				Collections.sort(match);
				matches.add(match);
				match = null;
			}
		}
		return matches;
	}

	private void loadClusters(int sequenceIdentity) {
		// load clusters only once
		if (clusters.size() > 0) {
			return;
		}

		if (!seqIdentities.contains(sequenceIdentity)) {
			logger.error("Representative chains are not available for %sequence identity: {}", sequenceIdentity);
			return;
		}

        String urlString = coreUrl + "bc-" + sequenceIdentity + ".out";

		try {

            URL u = new URL(urlString);
            InputStream stream = u.openStream();

            if (stream != null) {
                BufferedReader reader = new BufferedReader(new InputStreamReader(stream));

                String line = null;
                while ((line = reader.readLine()) != null) {
                    line = line.replaceAll("_", ".");
                    List<String> cluster = Arrays.asList(line.split(" "));
                    clusters.add(cluster);
                }
                reader.close();
                stream.close();
            } else {
                throw new IOException("Got null stream for URL " + urlString);
            }
        } catch (IOException e) {
		    logger.error("Could not get sequence clusters from URL " + urlString + ". Error: " + e.getMessage());
        }

	}

}

