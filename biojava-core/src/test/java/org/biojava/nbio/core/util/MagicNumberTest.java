package org.biojava.nbio.core.util;

import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

public class MagicNumberTest {

	@Test
	public void gzip() throws IOException {
		String file = this.getClass().getResource("example.gz").getFile();
		Path path = Paths.get(file);
		Assert.assertTrue("GZIP file", MagicNumber.isGZIP(path));

		file = this.getClass().getResource("build.xml").getFile();
		path = Paths.get(file);
		Assert.assertFalse("Not a GZIP file", MagicNumber.isGZIP(path));
	}
}
