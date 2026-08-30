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
 * Fetches the map coefficients published alongside the wwPDB validation reports.
 * <p>
 * <b>These are not density maps.</b> They are structure-factor amplitudes and
 * phases in mmCIF, exactly as used to produce the pictures in a validation
 * report, and a Fourier transform is required before anything can be drawn from
 * them &mdash; <code>gemmi sf2map</code>, or <code>cif2mtz</code> followed by
 * CCP4's <code>fft</code>. Nothing in BioJava or Jmol will render them.
 * <p>
 * They are supported because this is the route RCSB documents since
 * <code>edmaps.rcsb.org</code> was shut down in October 2024, and because they
 * are the authoritative archival form. For anything that needs to be displayed,
 * prefer {@link PdbeCcp4MapProvider} or {@link VolumeServerProvider}. Accordingly
 * this source is disabled by default in {@link DensityMapCache}.
 * <p>
 * One useful property: these servers return the content MD5 as the HTTP
 * <code>ETag</code>, so downloads from here are checksum-verified automatically.
 * <p>
 * The URLs are built against the documented download endpoint, which resolves an
 * entry by file name, rather than against the divided archive path. This is the
 * only provider here that ever had a choice &mdash; the density servers and the
 * EMDB archive are addressed by identifier already &mdash; and it matters because
 * the PDB moves to extended identifiers and a per-entry directory layout in July
 * 2027. A name survives that move; a constructed directory path does not. Mirrors
 * that publish directories instead of an endpoint are still reachable, through
 * {@link #DIVIDED_TWO_FO_FC_TEMPLATE} and {@link #ENTRIES_TWO_FO_FC_TEMPLATE}.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class WwpdbMapCoefficientsProvider extends AbstractDensityMapProvider {

	/** Default base URL, the wwPDB validation report download endpoint. */
	public static final String DEFAULT_SERVER_URL = "https://files.wwpdb.org/validation/download/";

	/** An RCSB mirror serving byte-identical files. */
	public static final String RCSB_MIRROR_URL = "https://files.rcsb.org/validation/download/";

	/**
	 * The wwPDB beta archive, which already holds the re-organised content that
	 * replaces the current archive on 21 July 2027 and serves the same endpoint.
	 * <p>
	 * Deliberately not the default, despite being the newer archive. The wwPDB
	 * describes this host as transitional: on the cutover date the beta archive
	 * replaces the main one, after which the beta URL is supported by redirection
	 * for three years. So it is the hostname that needs changing, twice, whereas
	 * {@link #DEFAULT_SERVER_URL} becomes the new archive and needs changing never.
	 * <p>
	 * Its value is in testing. Because this host is the post-2027 content today, a
	 * request against it checks the endpoint against the archive as it will be,
	 * rather than against the archive as it is.
	 */
	public static final String BETA_SERVER_URL = "https://files-beta.wwpdb.org/validation/download/";

	/**
	 * An EBI mirror serving byte-identical files.
	 * <p>
	 * Unlike the two above, EBI publishes no name-resolving endpoint &mdash; only
	 * full directory paths &mdash; so selecting it means setting the divided
	 * templates as well:
	 * <pre>
	 * setServerBaseUrl(EBI_MIRROR_URL);
	 * setPathUrlTemplate(DensityMapKind.TWO_FO_FC, DIVIDED_TWO_FO_FC_TEMPLATE);
	 * setPathUrlTemplate(DensityMapKind.FO_FC, DIVIDED_FO_FC_TEMPLATE);
	 * </pre>
	 * Setting the base alone yields 404s, because the flat file names do not exist
	 * there.
	 */
	public static final String EBI_MIRROR_URL = "https://ftp.ebi.ac.uk/pub/databases/pdb/validation_reports/";

	/**
	 * Default path template for the 2mFo-DFc coefficients: the file name alone.
	 * <p>
	 * The endpoint resolves an entry by name, so no directory path is built here.
	 * That is deliberate. In July 2027 the archive moves to extended identifiers
	 * and a per-entry directory layout, and a path assembled from a hash and an
	 * identifier would have to be rewritten for it; a name does not. Both spellings
	 * of an identifier resolve, so whichever {@code PdbId} yields is accepted.
	 */
	public static final String DEFAULT_TWO_FO_FC_TEMPLATE =
			"{pdbid_lc}_validation_2fo-fc_map_coef.cif.gz";

	/** Default path template for the mFo-DFc coefficients; see {@link #DEFAULT_TWO_FO_FC_TEMPLATE}. */
	public static final String DEFAULT_FO_FC_TEMPLATE =
			"{pdbid_lc}_validation_fo-fc_map_coef.cif.gz";

	/**
	 * Path template for the 2mFo-DFc coefficients in the divided archive, for
	 * mirrors that publish directories rather than an endpoint.
	 */
	public static final String DIVIDED_TWO_FO_FC_TEMPLATE =
			"{mid}/{pdbid_lc}/{pdbid_lc}_validation_2fo-fc_map_coef.cif.gz";

	/** Path template for the mFo-DFc coefficients in the divided archive. */
	public static final String DIVIDED_FO_FC_TEMPLATE =
			"{mid}/{pdbid_lc}/{pdbid_lc}_validation_fo-fc_map_coef.cif.gz";

	/**
	 * Path template for the 2mFo-DFc coefficients in the per-entry archive that
	 * replaces the divided one in July 2027, relative to a base URL ending in
	 * <code>.../pdb/data/</code>.
	 * <p>
	 * Provided so that a mirror of the new layout can be used the day it exists,
	 * without waiting for a release.
	 */
	public static final String ENTRIES_TWO_FO_FC_TEMPLATE =
			"entries/{mid}/{extid}/validation_reports/{extid}_validation_2fo-fc_map_coef.cif.gz";

	/** Path template for the mFo-DFc coefficients in the per-entry archive. */
	public static final String ENTRIES_FO_FC_TEMPLATE =
			"entries/{mid}/{extid}/validation_reports/{extid}_validation_fo-fc_map_coef.cif.gz";

	private static String serverBaseUrl = DEFAULT_SERVER_URL;
	private static String twoFoFcTemplate = DEFAULT_TWO_FO_FC_TEMPLATE;
	private static String foFcTemplate = DEFAULT_FO_FC_TEMPLATE;

	/**
	 * @param cacheRoot the BioJava cache directory
	 */
	public WwpdbMapCoefficientsProvider(File cacheRoot) {
		super(cacheRoot);
	}

	/** @return the base URL of the validation report archive */
	public static String getServerBaseUrl() {
		return serverBaseUrl;
	}

	/** @param url the base URL; a trailing slash is added if missing */
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
			throw new IllegalArgumentException("Map coefficients exist only for 2Fo-Fc and Fo-Fc, not " + kind);
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
		return DensityMapSource.WWPDB_MAP_COEFFICIENTS;
	}

	@Override
	public DensityFileFormat getFormat() {
		return DensityFileFormat.MAP_COEFFICIENTS_CIF_GZ;
	}

	@Override
	public boolean supports(DensityMapKind kind) {
		return kind == DensityMapKind.TWO_FO_FC || kind == DensityMapKind.FO_FC;
	}

	/**
	 * Builds the URL for a set of coefficients without fetching them.
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
