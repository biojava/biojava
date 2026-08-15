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
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

/**
 * Fetches downsampled volume slices from a Mol* density server.
 * <p>
 * This is usually the most economical source by a wide margin. For PDB entry
 * 1cbs the coarsest slice is about 210&nbsp;kB against roughly 2.1&nbsp;MB for
 * the two pre-computed CCP4 files, and for cryo-EM the difference is far larger
 * still: about 480&nbsp;kB against the 106&nbsp;MB primary map of EMD-0262. A
 * single response carries both the 2Fo-Fc and Fo-Fc data blocks for X-ray
 * entries.
 * <p>
 * <b>Caveats worth knowing.</b> Responses are chunked, with no
 * <code>Content-Length</code>, no <code>ETag</code> and no
 * <code>Last-Modified</code>, so the usual server-provided validation metadata is
 * unavailable; the size recorded for a cached slice is the byte count observed
 * during our own download instead. Detail levels are not linear either &mdash;
 * for a small entry the response stops growing past detail 1, while for a large
 * EM map the step from detail 2 to 3 multiplies the size several-fold.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class VolumeServerProvider extends AbstractDensityMapProvider {

	/**
	 * Which density server to talk to.
	 */
	public enum Host {
		/** RCSB's server at <code>maps.rcsb.org</code>. */
		RCSB(DensityMapSource.RCSB_VOLUME_SERVER, "https://maps.rcsb.org/", 3),
		/** PDBe's server at <code>www.ebi.ac.uk/pdbe/volume-server</code>. */
		PDBE(DensityMapSource.PDBE_VOLUME_SERVER, "https://www.ebi.ac.uk/pdbe/volume-server/", 6);

		private final DensityMapSource source;
		private final String defaultBaseUrl;
		private final int defaultDetail;

		Host(DensityMapSource source, String defaultBaseUrl, int defaultDetail) {
			this.source = source;
			this.defaultBaseUrl = defaultBaseUrl;
			this.defaultDetail = defaultDetail;
		}

		/** @return the source constant this host corresponds to */
		public DensityMapSource getSource() {
			return source;
		}

		/** @return the default base URL */
		public String getDefaultBaseUrl() {
			return defaultBaseUrl;
		}

		/** @return the detail level this host's own clients use by default */
		public int getDefaultDetail() {
			return defaultDetail;
		}
	}

	/** Default path template for an X-ray entry's whole unit cell. */
	public static final String DEFAULT_XRAY_CELL_TEMPLATE = "x-ray/{pdbid_lc}/cell?detail={detail}&encoding={encoding}";

	/** Default path template for an EM entry's whole cell. */
	public static final String DEFAULT_EM_CELL_TEMPLATE = "em/emd-{emdb_num}/cell?detail={detail}&encoding={encoding}";

	private final Host host;
	private String baseUrl;
	private int detail;
	private String encoding = "bcif";

	/**
	 * @param cacheRoot the BioJava cache directory
	 * @param host which density server to use
	 */
	public VolumeServerProvider(File cacheRoot, Host host) {
		super(cacheRoot);
		this.host = host;
		this.baseUrl = host.getDefaultBaseUrl();
		this.detail = host.getDefaultDetail();
	}

	/** @return which density server this instance uses */
	public Host getHost() {
		return host;
	}

	/** @return the base URL in use */
	public String getBaseUrl() {
		return baseUrl;
	}

	/** @param baseUrl the base URL; a trailing slash is added if missing */
	public void setBaseUrl(String baseUrl) {
		this.baseUrl = baseUrl == null ? host.getDefaultBaseUrl()
				: (baseUrl.endsWith("/") ? baseUrl : baseUrl + "/");
	}

	/** @return the detail level requested from the server */
	public int getDetail() {
		return detail;
	}

	/**
	 * Sets the detail level. Higher means a finer grid and a larger download; the
	 * relationship is neither linear nor the same for every entry.
	 *
	 * @param detail the level, normally 0 to 6
	 */
	public void setDetail(int detail) {
		this.detail = detail;
	}

	/** @return the encoding requested, either <code>bcif</code> or <code>cif</code> */
	public String getEncoding() {
		return encoding;
	}

	/**
	 * @param encoding <code>bcif</code> for BinaryCIF (default, and much smaller) or
	 *        <code>cif</code> for text
	 */
	public void setEncoding(String encoding) {
		this.encoding = encoding == null ? "bcif" : encoding;
	}

	@Override
	public DensityMapSource getSource() {
		return host.getSource();
	}

	@Override
	public DensityFileFormat getFormat() {
		return "cif".equalsIgnoreCase(encoding) ? DensityFileFormat.CIF_VOLUME : DensityFileFormat.BCIF_VOLUME;
	}

	@Override
	public boolean supports(DensityMapKind kind) {
		return kind == DensityMapKind.TWO_FO_FC || kind == DensityMapKind.FO_FC || kind == DensityMapKind.EM;
	}

	/**
	 * Builds the URL for a request without fetching it.
	 *
	 * @param request the request; its kind must be concrete
	 * @return the URL as a string, or <code>null</code> if the request cannot be served
	 */
	public String buildUrl(DensityMapRequest request) {
		if (request.getKind() == DensityMapKind.EM) {
			if (request.getEmdbId() == null) {
				return null;
			}
			return baseUrl + UrlTemplates.expand(withEncoding(DEFAULT_EM_CELL_TEMPLATE),
					UrlTemplates.values(null, request.getEmdbId(), detail));
		}
		if (request.getPdbId() == null) {
			return null;
		}
		return baseUrl + UrlTemplates.expand(withEncoding(DEFAULT_XRAY_CELL_TEMPLATE),
				UrlTemplates.values(urlId(request.getPdbId()), null, detail));
	}

	private String withEncoding(String template) {
		return template.replace("{encoding}", encoding);
	}

	@Override
	public DensityMapResult fetch(DensityMapRequest request) throws IOException {
		if (!supports(request.getKind())) {
			return null;
		}
		String urlString = buildUrl(request);
		if (urlString == null) {
			return null;
		}
		URL url = new URL(urlString);

		// The detail level changes the content, so it has to be part of the cache key.
		String qualifier = "d" + detail;
		File target;
		if (request.getKind() == DensityMapKind.EM) {
			target = DensityCacheLayout.emdbMapFile(effectiveCacheRoot(request), request.getEmdbId(),
					getSource(), getFormat(), qualifier);
		} else {
			// One response carries both the 2Fo-Fc and the Fo-Fc blocks, so the two
			// kinds share a cache entry: asking for each separately would otherwise
			// download and store the identical file twice. Which block to read is
			// decided at display time, not here.
			target = DensityCacheLayout.pdbMapFile(effectiveCacheRoot(request), request.getPdbId(),
					DensityCacheLayout.BOTH_KINDS_TOKEN, getSource(), getFormat(), qualifier);
		}
		DensityMapResult result = obtain(request, url, target, request.getKind(), request.getEmdbId(), null, null);

		if (request.getKind() == DensityMapKind.FO_FC) {
			return presentAsDifferenceMap(result);
		}
		return result;
	}

	/**
	 * Points the result at the companion file name that makes Jmol read the
	 * difference-map block. See
	 * {@link DensityCacheLayout#differenceMarkerFile(File)} for why the marker has
	 * to be in the name.
	 * <p>
	 * The companion is a hard link where the filesystem allows one, so the second
	 * name costs no additional space; a copy is only made if linking is refused.
	 */
	private DensityMapResult presentAsDifferenceMap(DensityMapResult result) throws IOException {
		File marker = DensityCacheLayout.differenceMarkerFile(result.getFile());
		if (!marker.isFile() || marker.length() != result.getFile().length()) {
			Files.deleteIfExists(marker.toPath());
			try {
				Files.createLink(marker.toPath(), result.getFile().toPath());
			} catch (IOException | UnsupportedOperationException e) {
				Files.copy(result.getFile().toPath(), marker.toPath(), StandardCopyOption.REPLACE_EXISTING);
			}
		}
		return new DensityMapResult(marker, result.getSource(), result.getFormat(), result.getKind(),
				result.getPdbId(), result.getEmdbId(), result.getSourceUrl(), result.isFromCache(),
				result.getRecommendedContourLevel(), result.getSigma());
	}
}
