package org.biojava.nbio.core.util;

import org.junit.Test;
import static org.junit.Assert.*;

import java.io.ByteArrayOutputStream;
import java.io.InputStream;

public class TestUncompressInputStream {
	
	/**
	 * The file compress_text.txt.lzc is the output of:
	 * <code>
	 *   cat compress_test.txt | compress > compress_test.txt.lzc
	 * </code>
	 * The original compress_test.txt contains text {@value #TEXT_IN_FILE} 
	 */
	private static final String TEST_FILE = "org/biojava/nbio/core/util/compress_test.txt.lzc";
	private static final String TEXT_IN_FILE = "Test of biojava uncompress.";
	
	@Test
	public void testUncompression() throws Exception {
		
		InputStream is = this.getClass().getClassLoader().getResourceAsStream(TEST_FILE);
		
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		UncompressInputStream.uncompress(is, baos);
		String decompressedText = baos.toString();
		
		assertEquals(TEXT_IN_FILE, decompressedText);
		
	}

}
