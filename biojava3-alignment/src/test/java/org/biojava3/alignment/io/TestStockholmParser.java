package org.biojava3.alignment.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

import org.biojava3.core.sequence.template.AbstractSequence;


import junit.framework.TestCase;

public class TestStockholmParser extends TestCase {

//	public void testStockholmParser(){
//
//		InputStream inStream = this.getClass().getResourceAsStream("/test.sth");
//
//		StockholmFileParser fileParser = new StockholmFileParser();
//
//		BufferedReader buf;
//		if (inStream == null) {
//			fail("input stream is null!");
//		}
//
//		buf = new BufferedReader(new InputStreamReader(inStream));
//
//		try {
//			StockholmStructure data = fileParser.parseFile(buf);					
//			System.out.println(data);
//
//			assertTrue(data.getBioSequences().size()==5);
//
//			AbstractSequence<?> seq = data.getBioSequences().get(0);
//			assertTrue(seq != null );
//
//			assertEquals(seq.getSequenceAsString(),"MTCRAQLIAVPRASSLAEAIACAQKMRVSRVPVYERS");
//
//		} catch (Exception e) {
//
//			e.printStackTrace();
//			fail(e.getMessage());
//		}
//
//
//	}

	public void testPiwi(){
		try {
			InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/piwi.sth.gz"));

			assertNotNull(inStream);

			StockholmFileParser fileParser = new StockholmFileParser();

			StockholmStructure data = fileParser.parse(inStream);					
			
			assertTrue("Did not get enough sequences!", data.getBioSequences(false).size()==20);

			AbstractSequence<?> seq = data.getBioSequences(false).get(0);
			assertTrue(seq != null );
			
			assertTrue(seq.getSequenceAsString().length() > 20);

		} catch (Exception e) {

			e.printStackTrace();
			fail(e.getMessage());
		}


	}

}
