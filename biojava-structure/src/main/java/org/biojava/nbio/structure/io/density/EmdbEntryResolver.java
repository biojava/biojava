/**
 * BioJava development code
 *
 * This code may be freely distributed and modified under the terms of the GNU
 * Lesser General Public Licence. This should be distributed with the code. If
 * you do not have a copy, see:
 *
 * http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual authors. These
 * should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims, or to join the
 * biojava-l mailing list, visit the home page at:
 *
 * http://www.biojava.org/
 */
package org.biojava.nbio.structure.io.density;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.time.Duration;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Properties;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;

import org.biojava.nbio.structure.PdbId;
import org.biojava.nbio.structure.align.util.URLConnectionTools;
import org.biojava.nbio.structure.io.LocalPDBDirectory.FetchBehavior;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Finds the EMDB entry, if any, associated with a PDB entry, and reads the few
 * facts about its map that are needed to fetch and display it.
 * <p>
 * The primary lookup is EMDB's own search API, which returns the identifier and
 * the author-recommended contour level together in a single small CSV response.
 * That is also the route Jmol uses, so BioJava's answer matches what a user sees
 * elsewhere. If it fails, RCSB's entry API is consulted for the identifier alone.
 * <p>
 * The experimental method is deliberately <i>not</i> inferred from a structure's
 * resolution: BioJava is known to mis-parse resolution for some cryo-EM entries
 * (biojava/biojava#1000), so the presence of an EMDB identifier is the reliable
 * signal.
 * <p>
 * Answers are cached on disk. Under
 * {@link FetchBehavior#LOCAL_ONLY} the cached answer is used whatever its age and
 * no connection is opened.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class EmdbEntryResolver {

	private static final Logger logger = LoggerFactory.getLogger(EmdbEntryResolver.class);

	/**
	 * Default template for the EMDB search that maps a PDB entry to the EM
	 * reconstructions it was fitted into, returning the contour level as well.
	 */
	public static final String DEFAULT_SEARCH_URL_TEMPLATE =
			"https://www.ebi.ac.uk/emdb/api/search/fitted_pdbs:{pdbid_lc}?fl=emdb_id,map_contour_level_value&wt=csv";

	/** Default template for the EMDB entry map metadata. */
	public static final String DEFAULT_MAP_INFO_URL_TEMPLATE =
			"https://www.ebi.ac.uk/emdb/api/entry/map/{emdb_id}";

	/** Default template for the RCSB entry record, used as a fallback. */
	public static final String DEFAULT_RCSB_ENTRY_URL_TEMPLATE =
			"https://data.rcsb.org/rest/v1/core/entry/{pdbid_lc}";

	private static final int TIMEOUT_MILLIS = 30000;
	private static final ObjectMapper MAPPER = new ObjectMapper();

	private static String searchUrlTemplate = DEFAULT_SEARCH_URL_TEMPLATE;
	private static String mapInfoUrlTemplate = DEFAULT_MAP_INFO_URL_TEMPLATE;
	private static String rcsbEntryUrlTemplate = DEFAULT_RCSB_ENTRY_URL_TEMPLATE;

	private File cacheRoot;
	private FetchBehavior fetchBehavior = FetchBehavior.FETCH_FILES;
	private int mappingMaxAgeDays = 30;

	/**
	 * @param cacheRoot the BioJava cache directory
	 */
	public EmdbEntryResolver(File cacheRoot) {
		this.cacheRoot = cacheRoot;
	}

	/** @param cacheRoot the BioJava cache directory */
	public void setCacheRoot(File cacheRoot) {
		this.cacheRoot = cacheRoot;
	}

	/** @param fetchBehavior how aggressively to refresh cached answers */
	public void setFetchBehavior(FetchBehavior fetchBehavior) {
		this.fetchBehavior = fetchBehavior == null ? FetchBehavior.FETCH_FILES : fetchBehavior;
	}

	/**
	 * @param days how long a cached mapping stays fresh. Mappings change only when
	 *        an entry is re-released, so this can be generous.
	 */
	public void setMappingMaxAgeDays(int days) {
		this.mappingMaxAgeDays = days;
	}

	/** @param template the EMDB search URL template */
	public static void setSearchUrlTemplate(String template) {
		searchUrlTemplate = template == null ? DEFAULT_SEARCH_URL_TEMPLATE : template;
	}

	/** @param template the EMDB map metadata URL template */
	public static void setMapInfoUrlTemplate(String template) {
		mapInfoUrlTemplate = template == null ? DEFAULT_MAP_INFO_URL_TEMPLATE : template;
	}

	/** @param template the RCSB entry URL template used as a fallback */
	public static void setRcsbEntryUrlTemplate(String template) {
		rcsbEntryUrlTemplate = template == null ? DEFAULT_RCSB_ENTRY_URL_TEMPLATE : template;
	}

	/** Restores all default URL templates. */
	public static void resetToDefaults() {
		searchUrlTemplate = DEFAULT_SEARCH_URL_TEMPLATE;
		mapInfoUrlTemplate = DEFAULT_MAP_INFO_URL_TEMPLATE;
		rcsbEntryUrlTemplate = DEFAULT_RCSB_ENTRY_URL_TEMPLATE;
	}

	/**
	 * The EMDB entries a PDB entry was fitted into.
	 *
	 * @param pdbId the PDB entry
	 * @return the EMDB identifiers, most relevant first; empty if the entry is not
	 *         an EM structure or has no associated map
	 */
	public List<String> getEmdbIds(PdbId pdbId) {
		Mapping cached = readMapping(pdbId);
		if (cached != null && (fetchBehavior == FetchBehavior.LOCAL_ONLY || cached.isFresh(mappingMaxAgeDays))) {
			return cached.emdbIds;
		}
		if (fetchBehavior == FetchBehavior.LOCAL_ONLY) {
			return Collections.emptyList();
		}

		Mapping fetched = searchEmdb(pdbId);
		if (fetched == null) {
			fetched = queryRcsb(pdbId);
		}
		if (fetched == null) {
			// Serve a stale answer rather than nothing when the services are unreachable.
			return cached == null ? Collections.emptyList() : cached.emdbIds;
		}
		writeMapping(pdbId, fetched);
		return fetched.emdbIds;
	}

	/**
	 * The contour level the depositors recommend for the map a PDB entry was fitted
	 * into, as reported by the EMDB search.
	 *
	 * @param pdbId the PDB entry
	 * @return the level in absolute map units, or <code>null</code> if unknown
	 */
	public Double getContourLevelFor(PdbId pdbId) {
		Mapping cached = readMapping(pdbId);
		if (cached != null && (fetchBehavior == FetchBehavior.LOCAL_ONLY || cached.isFresh(mappingMaxAgeDays))) {
			return cached.contourLevel;
		}
		getEmdbIds(pdbId);
		Mapping refreshed = readMapping(pdbId);
		return refreshed == null ? null : refreshed.contourLevel;
	}

	/**
	 * Reads an EMDB entry's map metadata: contour level, RMS deviation and the size
	 * of the primary map file.
	 * <p>
	 * The response is cached verbatim, so asking twice costs one request.
	 *
	 * @param emdbId the EMDB entry, in any accepted form
	 * @return the metadata, or <code>null</code> if it could not be obtained
	 */
	public EmdbEntryInfo getEntryInfo(String emdbId) {
		String canonical = DensityMapRequest.normalizeEmdbId(emdbId);
		File cacheFile = DensityCacheLayout.emdbMapInfoFile(cacheRoot, canonical);

		String json = null;
		if (cacheFile.isFile()) {
			try {
				json = new String(Files.readAllBytes(cacheFile.toPath()), StandardCharsets.UTF_8);
			} catch (IOException e) {
				logger.debug("Could not read cached EMDB metadata [{}]: {}", cacheFile, e.getMessage());
			}
		}
		if (json == null) {
			if (fetchBehavior == FetchBehavior.LOCAL_ONLY) {
				return null;
			}
			String url = UrlTemplates.expand(mapInfoUrlTemplate, UrlTemplates.values(null, canonical, -1));
			try {
				json = read(new URL(url));
			} catch (IOException e) {
				logger.warn("Could not read EMDB metadata for {}: {}", canonical, e.getMessage());
				return null;
			}
			writeCache(cacheFile, json);
		}

		try {
			JsonNode map = MAPPER.readTree(json).path("map");
			Double contour = null;
			JsonNode contours = map.path("contour_list").path("contour");
			for (JsonNode c : contours) {
				if (contour == null || c.path("primary").asBoolean(false)) {
					contour = c.path("level").isNumber() ? c.path("level").asDouble() : contour;
				}
			}
			Double sigma = map.path("statistics").path("std").isNumber()
					? map.path("statistics").path("std").asDouble() : null;
			Long bytes = map.path("size_kbytes").isNumber()
					? map.path("size_kbytes").asLong() * 1024L : null;
			return new EmdbEntryInfo(canonical, contour, sigma, bytes);
		} catch (IOException | RuntimeException e) {
			logger.warn("Could not parse EMDB metadata for {}: {}", canonical, e.getMessage());
			return null;
		}
	}

	private Mapping searchEmdb(PdbId pdbId) {
		String url = UrlTemplates.expand(searchUrlTemplate,
				UrlTemplates.values(DensityCacheLayout.shortIdOrFull(pdbId), null, -1));
		try {
			String csv = read(new URL(url));
			List<String> ids = new ArrayList<>();
			Double contour = null;
			String[] lines = csv.split("\\R");
			for (int i = 1; i < lines.length; i++) { // line 0 is the header
				String line = lines[i].trim();
				if (line.isEmpty()) {
					continue;
				}
				String[] cols = line.split(",");
				if (cols.length > 0 && !cols[0].isEmpty()) {
					ids.add(DensityMapRequest.normalizeEmdbId(cols[0]));
				}
				if (contour == null && cols.length > 1 && !cols[1].isEmpty()) {
					try {
						contour = Double.valueOf(cols[1].trim());
					} catch (NumberFormatException ignored) {
						// the column is optional and occasionally blank
					}
				}
			}
			return new Mapping(ids, contour, Instant.now());
		} catch (IOException e) {
			logger.debug("EMDB search for {} failed: {}", pdbId.getId(), e.getMessage());
			return null;
		}
	}

	private Mapping queryRcsb(PdbId pdbId) {
		String url = UrlTemplates.expand(rcsbEntryUrlTemplate,
				UrlTemplates.values(DensityCacheLayout.shortIdOrFull(pdbId), null, -1));
		try {
			JsonNode root = MAPPER.readTree(read(new URL(url)));
			JsonNode ids = root.path("rcsb_entry_container_identifiers").path("emdb_ids");
			List<String> result = new ArrayList<>();
			for (JsonNode id : ids) {
				result.add(DensityMapRequest.normalizeEmdbId(id.asText()));
			}
			return new Mapping(result, null, Instant.now());
		} catch (IOException | RuntimeException e) {
			logger.debug("RCSB entry lookup for {} failed: {}", pdbId.getId(), e.getMessage());
			return null;
		}
	}

	private String read(URL url) throws IOException {
		try (InputStream in = URLConnectionTools.getInputStream(url, true, TIMEOUT_MILLIS);
				BufferedReader reader = new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8))) {
			StringBuilder sb = new StringBuilder();
			char[] buffer = new char[8192];
			int n;
			while ((n = reader.read(buffer)) != -1) {
				sb.append(buffer, 0, n);
			}
			return sb.toString();
		}
	}

	private Mapping readMapping(PdbId pdbId) {
		File file = DensityCacheLayout.emdbMappingFile(cacheRoot, pdbId);
		if (!file.isFile()) {
			return null;
		}
		Properties p = new Properties();
		try (InputStream in = Files.newInputStream(file.toPath())) {
			p.load(in);
		} catch (IOException e) {
			return null;
		}
		String ids = p.getProperty("emdbIds", "");
		List<String> list = ids.isEmpty() ? Collections.emptyList() : Arrays.asList(ids.split(","));
		Double contour = null;
		try {
			String c = p.getProperty("contourLevel", "");
			contour = c.isEmpty() ? null : Double.valueOf(c);
		} catch (NumberFormatException ignored) {
			// leave it null
		}
		Instant retrieved;
		try {
			retrieved = Instant.parse(p.getProperty("retrieved"));
		} catch (RuntimeException e) {
			retrieved = Instant.EPOCH;
		}
		return new Mapping(list, contour, retrieved);
	}

	private void writeMapping(PdbId pdbId, Mapping mapping) {
		File file = DensityCacheLayout.emdbMappingFile(cacheRoot, pdbId);
		File dir = file.getParentFile();
		if (!dir.isDirectory() && !dir.mkdirs()) {
			logger.debug("Could not create [{}]", dir);
			return;
		}
		Properties p = new Properties();
		p.setProperty("pdbId", DensityCacheLayout.shortIdOrFull(pdbId));
		p.setProperty("emdbIds", String.join(",", mapping.emdbIds));
		p.setProperty("contourLevel", mapping.contourLevel == null ? "" : mapping.contourLevel.toString());
		p.setProperty("retrieved", mapping.retrieved.toString());
		try (OutputStream out = Files.newOutputStream(file.toPath())) {
			p.store(out, "BioJava PDB to EMDB mapping");
		} catch (IOException e) {
			logger.debug("Could not cache the EMDB mapping for {}: {}", pdbId.getId(), e.getMessage());
		}
	}

	private void writeCache(File file, String content) {
		File dir = file.getParentFile();
		if (!dir.isDirectory() && !dir.mkdirs()) {
			return;
		}
		try {
			Files.write(file.toPath(), content.getBytes(StandardCharsets.UTF_8));
		} catch (IOException e) {
			logger.debug("Could not cache [{}]: {}", file, e.getMessage());
		}
	}

	/** A cached PDB-to-EMDB answer. */
	private static final class Mapping {
		final List<String> emdbIds;
		final Double contourLevel;
		final Instant retrieved;

		Mapping(List<String> emdbIds, Double contourLevel, Instant retrieved) {
			this.emdbIds = Collections.unmodifiableList(new ArrayList<>(emdbIds));
			this.contourLevel = contourLevel;
			this.retrieved = retrieved;
		}

		boolean isFresh(int maxAgeDays) {
			return Duration.between(retrieved, Instant.now()).toDays() < maxAgeDays;
		}
	}
}
