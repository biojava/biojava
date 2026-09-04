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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;

import org.biojava.nbio.structure.PdbId;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Test;

/**
 * URL construction for each provider. These run offline: only the strings are
 * built, nothing is fetched.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class TestDensityMapUrlTemplates {

	private static final File ROOT = new File("/tmp/bjcache");

	@AfterEach
	public void restoreDefaults() {
		PdbeCcp4MapProvider.resetToDefaults();
		WwpdbMapCoefficientsProvider.resetToDefaults();
		EmdbMapProvider.resetToDefaults();
	}

	/**
	 * PDBe answers only to the lower-case four-character spelling; an upper-case or
	 * extended one returns HTTP 404.
	 */
	@Test
	public void pdbeUrlsAreLowerCase() {
		PdbeCcp4MapProvider p = new PdbeCcp4MapProvider(ROOT);
		assertEquals("https://www.ebi.ac.uk/pdbe/coordinates/files/1cbs.ccp4",
				p.buildUrl(new PdbId("1CBS"), DensityMapKind.TWO_FO_FC));
		assertEquals("https://www.ebi.ac.uk/pdbe/coordinates/files/1cbs_diff.ccp4",
				p.buildUrl(new PdbId("1cbs"), DensityMapKind.FO_FC));
	}

	@Test
	public void pdbeServerIsConfigurable() {
		PdbeCcp4MapProvider.setServerBaseUrl("http://localhost:1/maps");
		PdbeCcp4MapProvider p = new PdbeCcp4MapProvider(ROOT);
		assertEquals("http://localhost:1/maps/1cbs.ccp4",
				p.buildUrl(new PdbId("1cbs"), DensityMapKind.TWO_FO_FC));
	}

	/**
	 * The coefficients are addressed by file name, not by a constructed directory
	 * path, so that the July 2027 move to the per-entry archive does not invalidate
	 * the URLs.
	 */
	@Test
	public void wwpdbUrlsResolveByName() {
		WwpdbMapCoefficientsProvider p = new WwpdbMapCoefficientsProvider(ROOT);
		assertEquals("https://files.wwpdb.org/validation/download/"
				+ "1cbs_validation_2fo-fc_map_coef.cif.gz",
				p.buildUrl(new PdbId("1cbs"), DensityMapKind.TWO_FO_FC));
		assertEquals("https://files.wwpdb.org/validation/download/"
				+ "1cbs_validation_fo-fc_map_coef.cif.gz",
				p.buildUrl(new PdbId("1cbs"), DensityMapKind.FO_FC));
	}

	/**
	 * The identifier spelling is not hard-coded either way. An entry with a short
	 * form is asked for by it; one without &mdash; every entry deposited once the
	 * four-character space is exhausted &mdash; is asked for by its extended form.
	 * Both spellings resolve at the endpoint, so nothing here has to change in 2027.
	 */
	@Test
	public void wwpdbUrlsFollowWhicheverSpellingTheEntryHas() {
		WwpdbMapCoefficientsProvider p = new WwpdbMapCoefficientsProvider(ROOT);
		assertTrue(p.buildUrl(new PdbId("pdb_00001cbs"), DensityMapKind.TWO_FO_FC)
				.endsWith("/1cbs_validation_2fo-fc_map_coef.cif.gz"),
				"a shortable entry is asked for by its short name");
		assertTrue(p.buildUrl(new PdbId("pdb_00012abc"), DensityMapKind.TWO_FO_FC)
				.endsWith("/pdb_00012abc_validation_2fo-fc_map_coef.cif.gz"),
				"an entry with no short form is asked for by its extended name");
	}

	/**
	 * The beta host is a host swap and nothing more: the same endpoint, the same
	 * file name. That is what makes it a test target rather than a default - it
	 * holds the post-2027 content, so a request against it exercises the archive as
	 * it will be, while the default host is the one that survives the cutover.
	 */
	@Test
	public void betaHostChangesNothingButTheHost() {
		WwpdbMapCoefficientsProvider p = new WwpdbMapCoefficientsProvider(ROOT);
		String viaDefault = p.buildUrl(new PdbId("1cbs"), DensityMapKind.TWO_FO_FC);
		WwpdbMapCoefficientsProvider.setServerBaseUrl(WwpdbMapCoefficientsProvider.BETA_SERVER_URL);
		String viaBeta = p.buildUrl(new PdbId("1cbs"), DensityMapKind.TWO_FO_FC);
		assertEquals(viaDefault.substring(viaDefault.lastIndexOf('/')),
				viaBeta.substring(viaBeta.lastIndexOf('/')));
		assertTrue(viaBeta.startsWith("https://files-beta.wwpdb.org/validation/download/"));
	}

	/** EBI publishes directories rather than an endpoint, so it needs the divided templates. */
	@Test
	public void dividedTemplatesStillReachMirrors() {
		WwpdbMapCoefficientsProvider.setServerBaseUrl(WwpdbMapCoefficientsProvider.EBI_MIRROR_URL);
		WwpdbMapCoefficientsProvider.setPathUrlTemplate(DensityMapKind.TWO_FO_FC,
				WwpdbMapCoefficientsProvider.DIVIDED_TWO_FO_FC_TEMPLATE);
		WwpdbMapCoefficientsProvider p = new WwpdbMapCoefficientsProvider(ROOT);
		assertEquals("https://ftp.ebi.ac.uk/pub/databases/pdb/validation_reports/"
				+ "cb/1cbs/1cbs_validation_2fo-fc_map_coef.cif.gz",
				p.buildUrl(new PdbId("1cbs"), DensityMapKind.TWO_FO_FC));
	}

	/**
	 * The archive layout that arrives in July 2027 is expressible as a template,
	 * so a mirror of it needs configuration rather than a new release.
	 */
	@Test
	public void perEntryTemplateUsesExtendedIdentifiers() {
		WwpdbMapCoefficientsProvider.setServerBaseUrl("https://files-beta.wwpdb.org/pub/wwpdb/pdb/data/");
		WwpdbMapCoefficientsProvider.setPathUrlTemplate(DensityMapKind.TWO_FO_FC,
				WwpdbMapCoefficientsProvider.ENTRIES_TWO_FO_FC_TEMPLATE);
		WwpdbMapCoefficientsProvider p = new WwpdbMapCoefficientsProvider(ROOT);
		assertEquals("https://files-beta.wwpdb.org/pub/wwpdb/pdb/data/"
				+ "entries/cb/pdb_00001cbs/validation_reports/"
				+ "pdb_00001cbs_validation_2fo-fc_map_coef.cif.gz",
				p.buildUrl(new PdbId("1cbs"), DensityMapKind.TWO_FO_FC));
	}

	/** The extended spelling is available to any template, in the case the archive uses. */
	@Test
	public void extendedIdentifierPlaceholder() {
		assertEquals("pdb_00001cbs", UrlTemplates.values("1cbs", null, -1).get("extid"));
		assertEquals("pdb_00001cbs", UrlTemplates.values("pdb_00001cbs", null, -1).get("extid"));
		assertEquals("{extid}",
				UrlTemplates.expand("{extid}", UrlTemplates.values("not an id", null, -1)),
				"an unusable identifier leaves the placeholder alone rather than guessing");
	}

	/**
	 * The two-character directory is counted from the right hand end, so it is the
	 * same for both spellings. This is what lets the cache layout stay as it is
	 * across the 2027 transition, and it is the rule the wwPDB documents for the
	 * new archive.
	 */
	@Test
	public void middleHashIsTheSameForBothSpellings() {
		assertEquals("cb", UrlTemplates.values("1cbs", null, -1).get("mid"));
		assertEquals("cb", UrlTemplates.values("pdb_00001cbs", null, -1).get("mid"));
	}

	@Test
	public void volumeServerUrlsCarryDetailAndEncoding() {
		VolumeServerProvider rcsb = new VolumeServerProvider(ROOT, VolumeServerProvider.Host.RCSB);
		rcsb.setDetail(0);
		String url = rcsb.buildUrl(DensityMapRequest.builder(new PdbId("1cbs"))
				.kind(DensityMapKind.TWO_FO_FC).build());
		assertEquals("https://maps.rcsb.org/x-ray/1cbs/cell?detail=0&encoding=bcif", url);

		VolumeServerProvider pdbe = new VolumeServerProvider(ROOT, VolumeServerProvider.Host.PDBE);
		pdbe.setDetail(6);
		assertEquals("https://www.ebi.ac.uk/pdbe/volume-server/x-ray/1cbs/cell?detail=6&encoding=bcif",
				pdbe.buildUrl(DensityMapRequest.builder(new PdbId("1cbs"))
						.kind(DensityMapKind.TWO_FO_FC).build()));
	}

	@Test
	public void volumeServerEmUrlUsesTheEmdbNumber() {
		VolumeServerProvider rcsb = new VolumeServerProvider(ROOT, VolumeServerProvider.Host.RCSB);
		rcsb.setDetail(3);
		String url = rcsb.buildUrl(DensityMapRequest.builder("EMD-0262").build());
		assertEquals("https://maps.rcsb.org/em/emd-0262/cell?detail=3&encoding=bcif", url);
	}

	@Test
	public void encodingSelectsTheFormat() {
		VolumeServerProvider p = new VolumeServerProvider(ROOT, VolumeServerProvider.Host.RCSB);
		assertEquals(DensityFileFormat.BCIF_VOLUME, p.getFormat());
		p.setEncoding("cif");
		assertEquals(DensityFileFormat.CIF_VOLUME, p.getFormat());
		assertTrue(p.buildUrl(DensityMapRequest.builder(new PdbId("1cbs"))
				.kind(DensityMapKind.TWO_FO_FC).build()).endsWith("encoding=cif"));
	}

	@Test
	public void emdbMapUrl() {
		EmdbMapProvider p = new EmdbMapProvider(ROOT, null);
		assertEquals("https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0262/map/emd_0262.map.gz",
				p.buildUrl("EMD-0262"));
		assertEquals("https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0262/map/emd_0262.map.gz",
				p.buildUrl("emd_0262"));
	}

	/** An unknown placeholder must survive verbatim, so a bad template is obvious. */
	@Test
	public void unknownPlaceholdersAreLeftAlone() {
		assertEquals("a/{nosuch}/b",
				UrlTemplates.expand("a/{nosuch}/b", UrlTemplates.values("1cbs", null, -1)));
	}

	@Test
	public void providersDeclareWhatTheyCanServe() {
		assertTrue(new PdbeCcp4MapProvider(ROOT).supports(DensityMapKind.TWO_FO_FC));
		assertTrue(!new PdbeCcp4MapProvider(ROOT).supports(DensityMapKind.EM));
		assertTrue(new EmdbMapProvider(ROOT, null).supports(DensityMapKind.EM));
		assertTrue(!new EmdbMapProvider(ROOT, null).supports(DensityMapKind.TWO_FO_FC));
		assertTrue(new VolumeServerProvider(ROOT, VolumeServerProvider.Host.RCSB).supports(DensityMapKind.EM));
	}

	/** Only the coefficients are unrenderable; every sampled grid format is fine. */
	@Test
	public void onlyCoefficientsAreUnrenderable() {
		assertTrue(!DensityFileFormat.MAP_COEFFICIENTS_CIF_GZ.isJmolLoadable());
		assertTrue(DensityFileFormat.CCP4.isJmolLoadable());
		assertTrue(DensityFileFormat.CCP4_GZ.isJmolLoadable());
		assertTrue(DensityFileFormat.BCIF_VOLUME.isJmolLoadable());
		assertTrue(DensityFileFormat.CIF_VOLUME.isJmolLoadable());
	}
}
