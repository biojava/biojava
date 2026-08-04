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
package org.biojava.nbio.structure.test.util;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;

import org.biojava.nbio.structure.scop.ScopFactory;
import org.biojava.nbio.structure.scop.ScopInstallation;

/**
 * Installs small mock SCOP classification files from test resources so unit tests
 * do not depend on downloading from scop.berkeley.edu.
 */
public final class ScopTestHelper {

	private static final String RESOURCE_PREFIX = "scop/" + ScopInstallation.claFileName;

	private ScopTestHelper() {
	}

	/**
	 * Creates a {@link ScopInstallation} backed by a temp cache containing the mock
	 * {@code dir.cla.scop.txt_&lt;version&gt;} resource, and registers it with {@link ScopFactory}.
	 *
	 * @param version SCOP/SCOPe version, e.g. {@link ScopFactory#LATEST_VERSION}
	 * @return the installed database
	 */
	public static ScopInstallation installMockScop(String version) throws IOException {
		ScopInstallation scop = createMockScop(version);
		ScopFactory.setScopDatabase(scop);
		return scop;
	}

	/**
	 * Creates a {@link ScopInstallation} that reads a bundled mock classification file
	 * instead of downloading from the network.
	 *
	 * @param version SCOP/SCOPe version, e.g. {@link ScopFactory#VERSION_1_71}
	 * @return a SCOP database pointing at a temp cache with the mock file
	 */
	public static ScopInstallation createMockScop(String version) throws IOException {
		String resourceName = RESOURCE_PREFIX + version;
		InputStream in = ScopTestHelper.class.getClassLoader().getResourceAsStream(resourceName);
		if (in == null) {
			throw new IOException("Missing mock SCOP resource on classpath: " + resourceName);
		}

		Path cacheDir = Files.createTempDirectory("biojava-scop-mock-");
		cacheDir.toFile().deleteOnExit();

		Path claFile = cacheDir.resolve(ScopInstallation.claFileName + version);
		try (InputStream stream = in) {
			Files.copy(stream, claFile, StandardCopyOption.REPLACE_EXISTING);
		}
		claFile.toFile().deleteOnExit();

		ScopInstallation scop = new ScopInstallation(cacheDir.toAbsolutePath().toString());
		scop.setScopVersion(version);
		return scop;
	}
}
