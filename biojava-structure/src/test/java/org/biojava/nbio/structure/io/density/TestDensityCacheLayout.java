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
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.util.HashSet;
import java.util.Set;

import org.biojava.nbio.structure.PdbId;
import org.biojava.nbio.structure.io.LocalPDBDirectory;
import org.junit.jupiter.api.Test;

/**
 * Cache layout: directory derivation, and that no two source and kind
 * combinations can land on the same file.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class TestDensityCacheLayout {

	private static final File ROOT = new File("/tmp/bjcache");

	/**
	 * The divided-archive directory is taken from the end of the identifier, not
	 * the start. That is what makes both spellings of an entry agree: counting from
	 * the start would file the extended form under "db" instead of "cb".
	 */
	@Test
	public void middleHashIsTakenFromTheEndOfTheIdentifier() {
		assertEquals("cb", LocalPDBDirectory.getMiddleHash("1cbs"));
		assertEquals("cb", LocalPDBDirectory.getMiddleHash("1CBS"));
		assertEquals("cb", LocalPDBDirectory.getMiddleHash("pdb_00001cbs"));
		assertEquals("cb", LocalPDBDirectory.getMiddleHash("PDB_00001CBS"));

		// Guard against a regression to the substring(1, 3) form.
		assertFalse("db".equals(LocalPDBDirectory.getMiddleHash("pdb_00001cbs")));
	}

	@Test
	public void bothSpellingsOfAnEntryShareADirectory() {
		File shortForm = DensityCacheLayout.pdbMapFile(ROOT, new PdbId("1cbs"),
				DensityMapKind.TWO_FO_FC, DensityMapSource.PDBE_CCP4, DensityFileFormat.CCP4, null);
		File extendedForm = DensityCacheLayout.pdbMapFile(ROOT, new PdbId("PDB_00001CBS"),
				DensityMapKind.TWO_FO_FC, DensityMapSource.PDBE_CCP4, DensityFileFormat.CCP4, null);
		assertEquals(shortForm, extendedForm);
	}

	@Test
	public void pdbMapFileNamesAreLowerCaseAndDivided() {
		File f = DensityCacheLayout.pdbMapFile(ROOT, new PdbId("1CBS"),
				DensityMapKind.TWO_FO_FC, DensityMapSource.PDBE_CCP4, DensityFileFormat.CCP4, null);
		assertEquals("1cbs_2fofc_pdbe.ccp4", f.getName());
		assertEquals("cb", f.getParentFile().getName());
		assertEquals(DensityCacheLayout.DENSITY_DIR, f.getParentFile().getParentFile().getName());
	}

	@Test
	public void qualifierDistinguishesDetailLevels() {
		File d0 = DensityCacheLayout.pdbMapFile(ROOT, new PdbId("1cbs"), DensityMapKind.TWO_FO_FC,
				DensityMapSource.RCSB_VOLUME_SERVER, DensityFileFormat.BCIF_VOLUME, "d0");
		File d3 = DensityCacheLayout.pdbMapFile(ROOT, new PdbId("1cbs"), DensityMapKind.TWO_FO_FC,
				DensityMapSource.RCSB_VOLUME_SERVER, DensityFileFormat.BCIF_VOLUME, "d3");
		assertFalse(d0.equals(d3));
		assertTrue(d0.getName().endsWith("_d0.bcif"));
	}

	/**
	 * Every combination has to map to a distinct file, or one source would serve a
	 * cache hit belonging to another.
	 */
	@Test
	public void everySourceAndKindCombinationIsDistinct() {
		Set<String> seen = new HashSet<>();
		PdbId id = new PdbId("1cbs");
		for (DensityMapSource source : DensityMapSource.values()) {
			for (DensityMapKind kind : DensityMapKind.values()) {
				if (kind == DensityMapKind.AUTO || kind == DensityMapKind.EM) {
					continue; // AUTO is never cached; EM is keyed by EMDB id instead
				}
				for (DensityFileFormat format : DensityFileFormat.values()) {
					File f = DensityCacheLayout.pdbMapFile(ROOT, id, kind, source, format, null);
					assertTrue(seen.add(f.getPath()), "duplicate cache path: " + f);
				}
			}
		}
	}

	@Test
	public void emdbMapsAreKeyedByEmdbIdentifier() {
		File f = DensityCacheLayout.emdbMapFile(ROOT, "emd-262", DensityMapSource.EMDB_MAP,
				DensityFileFormat.CCP4_GZ, null);
		assertEquals("emd_262.map.gz", f.getName());
		assertEquals("EMD-262", f.getParentFile().getName());
		assertEquals(DensityCacheLayout.EMDB_DIR, f.getParentFile().getParentFile().getName());
	}

	@Test
	public void emdbIdentifiersAreNormalised() {
		assertEquals("EMD-0262", DensityMapRequest.normalizeEmdbId("EMD-0262"));
		assertEquals("EMD-0262", DensityMapRequest.normalizeEmdbId("emd-0262"));
		assertEquals("EMD-0262", DensityMapRequest.normalizeEmdbId("EMD_0262"));
		assertEquals("EMD-0262", DensityMapRequest.normalizeEmdbId("0262"));
		assertEquals("0262", DensityMapRequest.emdbNumber("EMD-0262"));
	}

	/**
	 * The marker has to sit inside the file name, before the extension. Jmol reads
	 * the difference block only when it finds that text in the name itself.
	 */
	@Test
	public void differenceMarkerGoesInsideTheName() {
		File plain = new File("/tmp/bjcache/density/cb/1cbs_both_rcsbvs_d3.bcif");
		File marked = DensityCacheLayout.differenceMarkerFile(plain);
		assertEquals("1cbs_both_rcsbvs_d3&diff=1.bcif", marked.getName());
		assertEquals(plain.getAbsoluteFile().getParentFile(), marked.getParentFile());
	}

	@Test
	public void emdbMappingFileIsDivided() {
		File f = DensityCacheLayout.emdbMappingFile(ROOT, new PdbId("6hu9"));
		assertEquals("6hu9.emdb.properties", f.getName());
		assertEquals("hu", f.getParentFile().getName());
		assertEquals(DensityCacheLayout.EMDB_MAPPING_DIR, f.getParentFile().getParentFile().getName());
	}
}
