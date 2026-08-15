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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;

import org.biojava.nbio.structure.PdbId;
import org.junit.After;
import org.junit.Test;

/**
 * URL construction for each provider. These run offline: only the strings are
 * built, nothing is fetched.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class TestDensityMapUrlTemplates {

	private static final File ROOT = new File("/tmp/bjcache");

	@After
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

	/** The validation archive is divided by the same two characters as the rest of the PDB. */
	@Test
	public void wwpdbUrlsUseTheDividedLayout() {
		WwpdbMapCoefficientsProvider p = new WwpdbMapCoefficientsProvider(ROOT);
		assertEquals("https://files.wwpdb.org/pub/pdb/validation_reports/"
				+ "cb/1cbs/1cbs_validation_2fo-fc_map_coef.cif.gz",
				p.buildUrl(new PdbId("1cbs"), DensityMapKind.TWO_FO_FC));
		assertEquals("https://files.wwpdb.org/pub/pdb/validation_reports/"
				+ "cb/1cbs/1cbs_validation_fo-fc_map_coef.cif.gz",
				p.buildUrl(new PdbId("1cbs"), DensityMapKind.FO_FC));
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
