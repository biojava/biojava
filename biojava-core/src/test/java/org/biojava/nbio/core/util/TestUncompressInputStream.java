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
package org.biojava.nbio.core.util;

import java.io.BufferedInputStream;
import org.junit.Test;
import static org.junit.Assert.*;

import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import org.junit.Assert;

public class TestUncompressInputStream {
	
	/**
	 * The file compress_text.txt.lzc is the output of:
	 * <code>
	 *   cat compress_test.txt | compress > compress_test.txt.lzc
	 * </code>
	 * The original compress_test.txt contains text {@value #TEXT_IN_FILE} 
	 */
	private static final String TEST_FILE = "org/biojava/nbio/core/util/compress_test.txt.Z";
	private static final String TEXT_IN_FILE = "Test of biojava uncompress.\n";
	
	private static final String BIGGER_TEST_FILE = "org/biojava/nbio/core/util/build-copy.xml.Z";
	private static final String ORIG_OF_BIGGER_TEST_FILE = "org/biojava/nbio/core/util/build.xml";
        
	@Test
	public void testUncompression() throws Exception {
		
		InputStream is = this.getClass().getClassLoader().getResourceAsStream(TEST_FILE);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		UncompressInputStream.uncompress(is, baos);
		String decompressedText = baos.toString();
		
		assertEquals(TEXT_IN_FILE, decompressedText);

		is = this.getClass().getClassLoader().getResourceAsStream(BIGGER_TEST_FILE);
		baos = new ByteArrayOutputStream();
		UncompressInputStream.uncompress(is, baos);

		ByteArrayOutputStream obaos = new ByteArrayOutputStream();
		try (BufferedInputStream oin = new BufferedInputStream(
				this.getClass().getClassLoader()
				.getResourceAsStream(ORIG_OF_BIGGER_TEST_FILE));) {
			byte[] buf = new byte[100000];
			int len;
			while ((len = oin.read(buf)) >= 0)
				obaos.write(buf, 0, len);
		}

		Assert.assertArrayEquals(baos.toByteArray(), obaos.toByteArray());
	}
}