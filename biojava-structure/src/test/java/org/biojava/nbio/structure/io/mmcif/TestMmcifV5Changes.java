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

import static org.junit.Assert.*;

import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Locale;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.junit.Before;
import org.junit.Test;

/**
 * Test date related changes in mmCIF 5.0 format.
 * 
 * @author Peter Rose
 * @author Aleix Lafita 
 * 
 */
public class TestMmcifV5Changes {
	
	private Structure s;
	
	@Before
	public void setup() throws IOException, StructureException {
		
		ClassLoader classLoader = this.getClass().getClassLoader();
		String file = classLoader.getResource("org/biojava/nbio/structure/io/mmcif/1stp_v5.cif").getPath();
		s = StructureIO.getStructure(file);
		
	}
	
    @Test
	public void testReleaseDate() throws ParseException {
    	
	    SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd",Locale.US);
        Date releaseDate = dateFormat.parse("1992-10-15");
	    assertEquals(releaseDate, s.getPDBHeader().getRelDate());
	}
    
	@Test
	public void testDepositionDate() throws ParseException {
		
	    SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd",Locale.US);
        Date depositionDate = dateFormat.parse("1992-03-12");
	    assertEquals(depositionDate, s.getPDBHeader().getDepDate());
	    
	}
	
	@Test
	public void testRevisionDate() throws ParseException {
		
	    SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd",Locale.US);
        Date depositionDate = dateFormat.parse("2011-07-13");
	    assertEquals(depositionDate, s.getPDBHeader().getModDate());
	    
	}

}
