package org.biojava.nbio.structure.io.mmcif;

import static org.junit.Assert.*;

import java.io.IOException;
import java.io.InputStream;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Locale;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.MMCIFFileReader;
import org.junit.Test;

public class TestMmcifV5Changes {

	/**
	 * Test date related changes in mmCIF 5.0 format
	 * @throws IOException
	 * @throws ParseException 
	 */
	
    @Test
	public void testReleaseDate() throws IOException, ParseException {
		Structure s = getStructure("/1stp_v50.cif");
	    SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd",Locale.US);
        Date releaseDate = dateFormat.parse("1992-10-15");
        // TODO uncomment the following line once getRelDate has been implemented
	    // assertEquals(releaseDate, s.getPDBHeader().getRelDate());    
	}
    
	@Test
	public void testDepositionDate() throws IOException, ParseException {
		Structure s = getStructure("/1stp_v50.cif");
	    SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd",Locale.US);
        Date depositionDate = dateFormat.parse("1992-03-12");
	    assertEquals(depositionDate, s.getPDBHeader().getDepDate());
	    
	}
	
	@Test
	public void testRevisionDate() throws IOException, ParseException {
		Structure s = getStructure("/1stp_v50.cif");
	    SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd",Locale.US);
        Date depositionDate = dateFormat.parse("2011-07-13");
	    assertEquals(depositionDate, s.getPDBHeader().getModDate());
	    
	}

	private Structure getStructure(String fileName) throws IOException{
		InputStream inStream = this.getClass().getResourceAsStream(fileName);
		assertNotNull(inStream);

		MMCIFFileReader reader = new MMCIFFileReader();

		return reader.getStructure(inStream) ;
	}
}
