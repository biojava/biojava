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
package org.biojava.nbio.structure.io.mmcif;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Locale;

import org.biojava.nbio.structure.PDBHeader;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

/**
 * Test parsing header information from MmCif files.
 * 
 * @author Anthony Bradley
 * @author Aleix Lafita
 *
 */
public class TestParseMmcifHeader {

	/**
	 * Test we can parse R-work and R-free effectively.
	 */
	@Test
	public void testRfactors() throws IOException, StructureException {

		AtomCache atomCache = new AtomCache();
		atomCache.setUseMmCif(true);

		Structure structure = atomCache.getStructure("4cup");

		PDBHeader pdbHeader = structure.getPDBHeader();
		// Check they are the same
		assertEquals(pdbHeader.getRfree(), 0.2078f, 0.000001f);
		assertEquals(pdbHeader.getRwork(), 0.1763f, 0.000001f);

	}

	/**
	 * Test parsing dates from MMCIF file version 4.
	 */
	@Test
	public void testDatesV4() throws IOException, StructureException, ParseException {
		
		ClassLoader classLoader = this.getClass().getClassLoader();
		String file4 = classLoader.getResource("org/biojava/nbio/structure/io/mmcif/1stp_v4.cif").getPath();
		Structure s = StructureIO.getStructure(file4);
		
		SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd",Locale.US);
		
		Date modDate = dateFormat.parse("2011-07-13");
		assertEquals(modDate, s.getPDBHeader().getModDate());

		Date releaseDate = dateFormat.parse("1992-10-15");
		assertEquals(releaseDate, s.getPDBHeader().getRelDate());
		
		Date depositionDate = dateFormat.parse("1992-03-12");
		assertEquals(depositionDate, s.getPDBHeader().getDepDate());
		
	}
	
	/**
	 * Test parsing dates from MMCIF file version 5.
	 */
	@Test
	public void testDatesV5() throws IOException, StructureException, ParseException {
		
		ClassLoader classLoader = this.getClass().getClassLoader();
		String file5 = classLoader.getResource("org/biojava/nbio/structure/io/mmcif/1stp_v5.cif").getPath();
		Structure s = StructureIO.getStructure(file5);
		
		SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd",Locale.US);

		Date modDate = dateFormat.parse("2011-07-13");
		assertEquals(modDate, s.getPDBHeader().getModDate());

		Date releaseDate = dateFormat.parse("1992-10-15");
		assertEquals(releaseDate, s.getPDBHeader().getRelDate());
		
		Date depositionDate = dateFormat.parse("1992-03-12");
		assertEquals(depositionDate, s.getPDBHeader().getDepDate());
		

	}

}
