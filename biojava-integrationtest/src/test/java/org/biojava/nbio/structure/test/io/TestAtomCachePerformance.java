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
package org.biojava.nbio.structure.test.io;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.LocalPDBDirectory.FetchBehavior;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

/**
 * A test to check the performance of AtomCache downloading
 *
 * By default it is ignored.
 * To execute use:
 * <pre>
 * mvn -Dtest=TestAtomCachePerformance test
 * </pre>
 *
 * @author duarte_j
 *
 */
public class TestAtomCachePerformance {

	private static final String[] PDB_IDS = {
		"1zjo",		"2dqc",		"4af2",		"1r52",		"4f3u",		"1f9v",		"3kuq",		"2yr4",		"3m4f",		"4j5p",
		"7ccp",		"4kro",		"1x7q",		"2gaw",		"2kli",		"2bdo",		"3csf",		"1muu",		"190l",		"2ecm"
	 	,
		"2f0y",		"3ind",		"3uu6",		"1p9j",		"1vm7",		"2y2c",		"2hez",		"1yrm",		"1yzx",		"1ps9",
		"3ue0",		"2o0o",		"2g59",		"4ees",		"2yfc",		"2anr",		"3cxk",		"2e7t",		"3kmh",		"3h00",
		"3gdm",		"1c0t",		"1fi0",		"2kqt",		"1ky8",		"169l",		"1z6h",		"1wbm",		"4g1j",		"1v3c",
		"2chm",		"4f0n",		"2vxb",		"2w0q",		"1g1n",		"3o6g",		"4eug",		"3nrm",		"3heo",		"4ewe",
		"2xjb",		"1vgj",		"3tpp",		"2gnl",		"3jpz",		"2pgt",		"1fn2",		"2h13",		"1xyj",		"1ds7",
		"2x93",		"4j5y",		"2bk2",		"1v83",		"4lj9",		"4ahc",		"1m34",		"1jo4",		"3flb",		"2cb2",
		"4k3p",		"1yq8",		"2h7z",		"2lbp",		"3vas",		"4jwn",		"2e47",		"3r43",		"3edd",		"3kss",
		"2dnk",		"1kg2",		"2pwh",		"1sjh",		"4cc0",		"3a7c",		"1o5a",		"4fu7",		"3hc4",		"3hoz"
		};

	private static AtomCache cache;

	@BeforeClass
	public static void setUpBeforeClass() {
		cache = new AtomCache();
		cache.setFetchBehavior(FetchBehavior.FORCE_DOWNLOAD);
	}

	@Ignore
	@Test
	public void testDownload() throws IOException, StructureException {
		System.out.println("Starting performance test for "+PDB_IDS.length+" PDB ids");
		long start = System.currentTimeMillis();
		for (String pdbId:PDB_IDS) {
			Structure cifS = getCifStructure(pdbId);
			Structure pdbS = getPdbStructure(pdbId);
			assertNotNull(cifS);
			assertNotNull(pdbS);
			assertEquals(pdbId,cifS.getPDBCode().toLowerCase());
			assertEquals(cifS.getPDBCode(),pdbS.getPDBCode());

			//System.out.print(".");

		}

		System.out.println();

		long end = System.currentTimeMillis();

		System.out.printf("Done in %5.1f s\n",(end-start)/1000.0);
	}

	private Structure getCifStructure(String pdbId) throws IOException, StructureException {
		cache.setUseMmCif(true);

		return cache.getStructure(pdbId);

	}

	private Structure getPdbStructure(String pdbId) throws IOException, StructureException {
		cache.setUseMmCif(false);

		return cache.getStructure(pdbId);

	}
}
