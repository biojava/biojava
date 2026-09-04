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

import org.biojava.nbio.structure.PdbId;

/**
 * Fetches PDBe's pre-computed CCP4 maps.
 * <p>
 * These are full-resolution sampled grids, one file per map kind, and are the
 * source Jmol itself uses for its built-in map shortcuts. They are larger than a
 * density-server slice &mdash; about 1&nbsp;MB per map for a small entry such as
 * 1cbs &mdash; but need no interpretation beyond contouring.
 * <p>
 * The service accepts only the lower-case four-character spelling of an
 * identifier; <code>1CBS.ccp4</code> and <code>pdb_00001cbs.ccp4</code> both
 * return HTTP 404. Entries deposited without structure factors, and cryo-EM
 * entries, have no maps here at all.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class PdbeCcp4MapProvider extends AbstractDensityMapProvider {

	/** Default base URL of the PDBe map service. */
	public static final String DEFAULT_SERVER_URL = "https://www.ebi.ac.uk/pdbe/coordinates/files/";

	/**
	 * An equivalent base URL serving the same files, kept in the documentation as a
	 * ready alternative should the primary one change.
	 */
	public static final String ALTERNATIVE_SERVER_URL = "https://www.ebi.ac.uk/pdbe/entry-files/";

	/** Default path template for the 2mFo-DFc map. */
	public static final String DEFAULT_TWO_FO_FC_TEMPLATE = "{pdbid_lc}.ccp4";

	/** Default path template for the mFo-DFc difference map. */
	public static final String DEFAULT_FO_FC_TEMPLATE = "{pdbid_lc}_diff.ccp4";

	private static String serverBaseUrl = DEFAULT_SERVER_URL;
	private static String twoFoFcTemplate = DEFAULT_TWO_FO_FC_TEMPLATE;
	private static String foFcTemplate = DEFAULT_FO_FC_TEMPLATE;

	/**
	 * @param cacheRoot the BioJava cache directory
	 */
	public PdbeCcp4MapProvider(File cacheRoot) {
		super(cacheRoot);
	}

	/** @return the base URL of the map service */
	public static String getServerBaseUrl() {
		return serverBaseUrl;
	}

	/**
	 * @param url the base URL of the map service; a trailing slash is added if missing
	 */
	public static void setServerBaseUrl(String url) {
		serverBaseUrl = url == null ? DEFAULT_SERVER_URL : (url.endsWith("/") ? url : url + "/");
	}

	/**
	 * Overrides the path template for a map kind.
	 *
	 * @param kind {@link DensityMapKind#TWO_FO_FC} or {@link DensityMapKind#FO_FC}
	 * @param template a template understood by {@link UrlTemplates}
	 */
	public static void setPathUrlTemplate(DensityMapKind kind, String template) {
		if (kind == DensityMapKind.TWO_FO_FC) {
			twoFoFcTemplate = template == null ? DEFAULT_TWO_FO_FC_TEMPLATE : template;
		} else if (kind == DensityMapKind.FO_FC) {
			foFcTemplate = template == null ? DEFAULT_FO_FC_TEMPLATE : template;
		} else {
			throw new IllegalArgumentException("PDBe CCP4 maps exist only for 2Fo-Fc and Fo-Fc, not " + kind);
		}
	}

	/** Restores the default server and templates. */
	public static void resetToDefaults() {
		serverBaseUrl = DEFAULT_SERVER_URL;
		twoFoFcTemplate = DEFAULT_TWO_FO_FC_TEMPLATE;
		foFcTemplate = DEFAULT_FO_FC_TEMPLATE;
	}

	@Override
	public DensityMapSource getSource() {
		return DensityMapSource.PDBE_CCP4;
	}

	@Override
	public DensityFileFormat getFormat() {
		return DensityFileFormat.CCP4;
	}

	@Override
	public boolean supports(DensityMapKind kind) {
		return kind == DensityMapKind.TWO_FO_FC || kind == DensityMapKind.FO_FC;
	}

	/**
	 * Builds the URL for a map without fetching it.
	 *
	 * @param pdbId the entry
	 * @param kind the kind of map
	 * @return the URL as a string
	 */
	public String buildUrl(PdbId pdbId, DensityMapKind kind) {
		String template = kind == DensityMapKind.FO_FC ? foFcTemplate : twoFoFcTemplate;
		return serverBaseUrl + UrlTemplates.expand(template, UrlTemplates.values(urlId(pdbId), null, -1));
	}

	@Override
	public DensityMapResult fetch(DensityMapRequest request) throws IOException {
		if (request.getPdbId() == null || !supports(request.getKind())) {
			return null;
		}
		URL url = new URL(buildUrl(request.getPdbId(), request.getKind()));
		File target = DensityCacheLayout.pdbMapFile(effectiveCacheRoot(request), request.getPdbId(),
				request.getKind(), getSource(), getFormat(), null);
		return obtain(request, url, target, request.getKind(), null, null, null);
	}
}
