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

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.zip.GZIPOutputStream;

import org.biojava.nbio.core.util.FileDownloadUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

/**
 * The CCP4 header check that keeps a server's error page out of the cache.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class TestCcp4Header {

	private File dir;

	@Before
	public void setUp() throws IOException {
		dir = Files.createTempDirectory("bj-ccp4").toFile();
	}

	@After
	public void tearDown() throws IOException {
		FileDownloadUtils.deleteDirectory(dir.toPath());
	}

	/** A minimal file carrying the stamp at the offset a real CCP4 header uses. */
	private static byte[] fakeMap() {
		byte[] bytes = new byte[2048];
		byte[] stamp = Ccp4Header.MAP_STAMP.getBytes(StandardCharsets.US_ASCII);
		System.arraycopy(stamp, 0, bytes, Ccp4Header.MAP_STAMP_OFFSET, stamp.length);
		return bytes;
	}

	private File write(String name, byte[] content) throws IOException {
		File f = new File(dir, name);
		Files.write(f.toPath(), content);
		return f;
	}

	@Test
	public void recognisesAMapByItsStamp() throws IOException {
		assertTrue(Ccp4Header.isCcp4(write("good.ccp4", fakeMap())));
	}

	@Test
	public void recognisesAGzippedMap() throws IOException {
		ByteArrayOutputStream buffer = new ByteArrayOutputStream();
		try (GZIPOutputStream gz = new GZIPOutputStream(buffer)) {
			gz.write(fakeMap());
		}
		assertTrue("EMDB serves its maps gzipped", Ccp4Header.isCcp4(write("good.map.gz", buffer.toByteArray())));
	}

	/**
	 * The case this check exists for: a server answering with an error page and an
	 * HTTP 200, which nothing else would catch.
	 */
	@Test
	public void rejectsAnHtmlErrorPage() throws IOException {
		StringBuilder html = new StringBuilder("<html><head><title>404 Not Found</title></head><body>");
		while (html.length() < 1500) {
			html.append("<p>The requested resource was not found on this server.</p>");
		}
		html.append("</body></html>");
		assertFalse(Ccp4Header.isCcp4(write("error.ccp4", html.toString().getBytes(StandardCharsets.UTF_8))));
	}

	@Test
	public void rejectsRandomBytesAndShortFiles() throws IOException {
		byte[] noise = new byte[2048];
		for (int i = 0; i < noise.length; i++) {
			noise[i] = (byte) (i * 31);
		}
		assertFalse(Ccp4Header.isCcp4(write("noise.ccp4", noise)));
		assertFalse("a file shorter than the header cannot be a map",
				Ccp4Header.isCcp4(write("tiny.ccp4", new byte[10])));
	}

	@Test
	public void quietVariantSwallowsUnreadableFiles() {
		assertFalse(Ccp4Header.isCcp4Quietly(new File(dir, "does-not-exist.ccp4")));
	}
}
