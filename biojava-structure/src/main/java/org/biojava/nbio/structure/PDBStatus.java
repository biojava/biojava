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
package org.biojava.nbio.structure;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.type.TypeFactory;
import org.biojava.nbio.structure.align.util.URLConnectionTools;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.*;

/**
 * Methods for getting the status of a PDB file (current, removed, unreleased)
 * and for accessing different versions of the structure.
 *
 * <p>
 * All methods query the
 * <a href="https://data.rcsb.org">
 * RCSB Data REST API</a>
 * <p>
 *
 * @author Spencer Bliven
 * @author Amr ALHOSSARY
 * @author Jose Duarte
 * @since 3.0.2
 */
public class PDBStatus {

	private static final Logger logger = LoggerFactory.getLogger(PDBStatus.class);

	public static final String DEFAULT_RCSB_DATA_API_SERVER = "data.rcsb.org";
	public static final String ALL_CURRENT_ENDPOINT = "https://%s/rest/v1/holdings/current/entry_ids";
	public static final String STATUS_ENDPOINT = "https://%s/rest/v1/holdings/status/%s";
	public static final String STATUS_LIST_ENDPOINT = "https://%s/rest/v1/holdings/status?ids=%s";

	/**
	 * Represents a simplified 3 state status of PDB IDs.
	 * @author Spencer Bliven
	 */
	public enum Status {
		// the simplified status enum in rcsb_repository_holdings_combined
		REMOVED,
		CURRENT,
		UNRELEASED;

		/**
		 * @throws IllegalArgumentException If the string is not recognized
		 */
		public static Status fromString(String statusStr) {
			if (statusStr == null) throw new IllegalArgumentException("Status string can't be null");
			if(statusStr.equalsIgnoreCase("REMOVED"))
				return Status.REMOVED;
			else if(statusStr.equalsIgnoreCase("CURRENT"))
				return Status.CURRENT;
			else if(statusStr.equalsIgnoreCase("UNRELEASED"))
				return Status.UNRELEASED;
			else {
				throw new IllegalArgumentException("Unable to parse status '"+statusStr+"'.");
			}
		}
	}

	/**
	 * Get the status of a PDB id.
	 *
	 * @param pdbId the id
	 * @return The status.
	 */
	public static Status getStatus(String pdbId) throws IOException {
		URL url = new URL(String.format(STATUS_ENDPOINT, DEFAULT_RCSB_DATA_API_SERVER, pdbId.toUpperCase()));
		ObjectMapper objectMapper = new ObjectMapper();
		JsonNode node = objectMapper.readValue(url.openStream(), JsonNode.class);
		return parseStatusRecord(node);
	}

	/**
	 * Get the status of a collection of PDB ids (in a single API query).
	 *
	 * @see #getStatus(String)
	 * @param pdbIds the ids
	 * @return The status array
	 */
	public static Status[] getStatus(String[] pdbIds) throws IOException {

		URL url = new URL(String.format(STATUS_LIST_ENDPOINT, DEFAULT_RCSB_DATA_API_SERVER, String.join(",", pdbIds)));

		List<Status> statuses = new ArrayList<>();

		ObjectMapper objectMapper = new ObjectMapper();
		JsonNode node = objectMapper.readValue(url.openStream(), JsonNode.class);

		if (node !=null && node.isArray()) {
			for (JsonNode record : node) {
				Status status = parseStatusRecord(record);
				statuses.add(status);
			}
		}

		if (statuses.size() != pdbIds.length) {
			logger.warn("RCSB status request was for {} ids, but {} were returned", pdbIds.length, statuses.size());
		}

		return statuses.toArray(new Status[0]);
	}

	private static Status parseStatusRecord(JsonNode jsonNode) {
		// e.g.
		// "rcsb_repository_holdings_combined": {
		//"id_code_replaced_by_latest": "4HHB",
		//"status": "REMOVED",
		//"status_code": "OBS"
		//},
		JsonNode rcsbRepoHoldingsNode = jsonNode.get("rcsb_repository_holdings_combined");
		return Status.fromString(rcsbRepoHoldingsNode.get("status").asText());
	}

	/**
	 * Gets the current version of a PDB ID.
	 *
	 * @param oldPdbId the id
	 * @return The replacement for oldPdbId, or null if none are found.
	 * If entry is current then the input PDB id is returned
	 */
	public static String getCurrent(String oldPdbId) throws IOException {
		URL url = new URL(String.format(STATUS_ENDPOINT, DEFAULT_RCSB_DATA_API_SERVER, oldPdbId.toUpperCase()));
		ObjectMapper objectMapper = new ObjectMapper();
		JsonNode node = objectMapper.readValue(url.openStream(), JsonNode.class);
		JsonNode rcsbRepoHoldingsNode = node.get("rcsb_repository_holdings_combined");
		Status st = Status.fromString(rcsbRepoHoldingsNode.get("status").asText());
		if (st == Status.REMOVED) {
			JsonNode replacedByNode = rcsbRepoHoldingsNode.get("id_code_replaced_by_latest");
			if (replacedByNode != null)
				return replacedByNode.asText();
			else
				return null;
		} else if (st == Status.CURRENT) {
			return oldPdbId;
		} else {
			return null;
		}

	}

	/**
	 * Returns all current PDB IDs
	 *
	 * @return a list of PDB IDs
	 * @throws IOException if a problem occurs retrieving the information
	 */
	public static SortedSet<String> getCurrentPDBIds() throws IOException {

		// Build REST query URL
		String urlStr = String.format(ALL_CURRENT_ENDPOINT, DEFAULT_RCSB_DATA_API_SERVER);
		URL u = new URL(urlStr);

		InputStream stream = URLConnectionTools.getInputStream(u, 60000);

		ObjectMapper objectMapper = new ObjectMapper();
		TypeFactory typeFactory = objectMapper.getTypeFactory();
		List<String> pdbIdList = objectMapper.readValue(stream, typeFactory.constructCollectionType(List.class, String.class));

		return new TreeSet<>(pdbIdList);
	}

	public static void main(String[] args) throws Exception {
		SortedSet<String> all = getCurrentPDBIds();
		System.out.println("Number of current PDB ids is: " + all.size());
	}
}
