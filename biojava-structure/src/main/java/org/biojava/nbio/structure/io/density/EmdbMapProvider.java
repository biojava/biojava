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

import java.io.File;
import java.io.IOException;
import java.net.URL;

/**
 * Fetches the full-resolution primary map of an EMDB entry.
 * <p>
 * These are gzipped CCP4/MRC files, and they are big: the primary map of
 * EMD-0262 is about 106&nbsp;MB, and gigabyte maps exist. For most display
 * purposes a slice from a density server is a far better trade &mdash; a few
 * hundred kilobytes for the same entry &mdash; which is why
 * {@link DensityMapCache} tries {@link VolumeServerProvider} first and only falls
 * back here. Use this source when the full sampling genuinely matters.
 * <p>
 * The size limit is checked against the size EMDB itself reports before any of
 * the body is transferred, so exceeding it costs one small metadata request
 * rather than a partial download.
 * <p>
 * Jmol recognises gzip from the file's magic bytes, so the cached
 * <code>.map.gz</code> can be handed to it without decompressing first.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class EmdbMapProvider extends AbstractDensityMapProvider {

	/** Default base URL of the EMDB archive. */
	public static final String DEFAULT_SERVER_URL = "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/";

	/** Default path template for an entry's primary map. */
	public static final String DEFAULT_MAP_TEMPLATE = "{emdb_id}/map/emd_{emdb_num}.map.gz";

	private static String serverBaseUrl = DEFAULT_SERVER_URL;
	private static String mapTemplate = DEFAULT_MAP_TEMPLATE;

	private final EmdbEntryResolver resolver;

	/**
	 * @param cacheRoot the BioJava cache directory
	 * @param resolver used to find the EMDB entry for a PDB entry and to read its
	 *        contour level and size
	 */
	public EmdbMapProvider(File cacheRoot, EmdbEntryResolver resolver) {
		super(cacheRoot);
		this.resolver = resolver;
	}

	/** @return the base URL of the EMDB archive */
	public static String getServerBaseUrl() {
		return serverBaseUrl;
	}

	/** @param url the base URL; a trailing slash is added if missing */
	public static void setServerBaseUrl(String url) {
		serverBaseUrl = url == null ? DEFAULT_SERVER_URL : (url.endsWith("/") ? url : url + "/");
	}

	/** @param template the path template for an entry's primary map */
	public static void setMapUrlTemplate(String template) {
		mapTemplate = template == null ? DEFAULT_MAP_TEMPLATE : template;
	}

	/** Restores the default server and template. */
	public static void resetToDefaults() {
		serverBaseUrl = DEFAULT_SERVER_URL;
		mapTemplate = DEFAULT_MAP_TEMPLATE;
	}

	@Override
	public DensityMapSource getSource() {
		return DensityMapSource.EMDB_MAP;
	}

	@Override
	public DensityFileFormat getFormat() {
		return DensityFileFormat.CCP4_GZ;
	}

	@Override
	public boolean supports(DensityMapKind kind) {
		return kind == DensityMapKind.EM;
	}

	/**
	 * Builds the URL of an entry's primary map without fetching it.
	 *
	 * @param emdbId the EMDB entry, in any accepted form
	 * @return the URL as a string
	 */
	public String buildUrl(String emdbId) {
		return serverBaseUrl + UrlTemplates.expand(mapTemplate, UrlTemplates.values(null, emdbId, -1));
	}

	@Override
	public DensityMapResult fetch(DensityMapRequest request) throws IOException {
		if (!supports(request.getKind())) {
			return null;
		}
		String emdbId = request.getEmdbId();
		if (emdbId == null) {
			return null;
		}

		EmdbEntryInfo info = resolver == null ? null : resolver.getEntryInfo(emdbId);

		// Decline before transferring anything: EMDB reports the map size in its
		// metadata, so an oversized map costs one small request rather than a
		// partial multi-hundred-megabyte download.
		long limit = effectiveMaxBytes(request);
		if (limit > 0 && info != null && info.getMapSizeBytes() != null && info.getMapSizeBytes() > limit) {
			URL url = new URL(buildUrl(emdbId));
			reportTooLarge(url, info.getMapSizeBytes(), limit, request);
			throw new DensityMapTooLargeException(url.toString(), info.getMapSizeBytes(), limit);
		}

		URL url = new URL(buildUrl(emdbId));
		File target = DensityCacheLayout.emdbMapFile(effectiveCacheRoot(request), emdbId,
				getSource(), getFormat(), null);
		return obtain(request, url, target, DensityMapKind.EM, emdbId,
				info == null ? null : info.getRecommendedContourLevel(),
				info == null ? null : info.getSigma());
	}
}
