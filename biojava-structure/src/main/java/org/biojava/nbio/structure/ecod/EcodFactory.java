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
package org.biojava.nbio.structure.ecod;

import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.cath.CathDatabase;
import org.biojava.nbio.structure.cath.CathInstallation;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Controls global {@link CathDatabase CathDatabases} being used.
 * Implements a multiton pattern through {@link #getCathDatabase(String)},
 * and a singleton pattern through {@link #getCathDatabase()}.
 * @author Spencer Bliven
 * @see ScopFactory
 * @see CathInstallation
 */
public class EcodFactory {

	private static Logger logger = LoggerFactory.getLogger(EcodFactory.class);

	public static String DEFAULT_VERSION = EcodInstallation.DEFAULT_VERSION;

	private static Map<String, EcodDatabase> versionedEcodDBs =
			Collections.synchronizedMap(new HashMap<String, EcodDatabase>());
	private static String defaultVersion = DEFAULT_VERSION;

	/**
	 * Returns the (singleton) database for the current default version
	 */
	public static EcodDatabase getEcodDatabase() {
		return getEcodDatabase(defaultVersion);
	}

	public static EcodDatabase getEcodDatabase(String version) {
		if( version == null )
			version = defaultVersion;

		synchronized(versionedEcodDBs) {
			EcodDatabase ecod = versionedEcodDBs.get(version.toLowerCase());
			if( ecod == null ) {
				logger.info("Creating new {}, version {}",EcodInstallation.class.getSimpleName(),version);
				String cacheDir = new UserConfiguration().getCacheFilePath();
				ecod = new EcodInstallation(cacheDir, version);
				versionedEcodDBs.put(version.toLowerCase(), ecod);

				// If the parsed version differed from that requested, add that too
				try {
					if( ! versionedEcodDBs.containsKey(ecod.getVersion().toLowerCase()) ) {
						versionedEcodDBs.put(ecod.getVersion().toLowerCase(),ecod);
					}
				} catch (IOException e) {
					// For parsing errors, just use the requested version
				}
			}

			return ecod;
		}
	}

	/**
	 * Updates the default version
	 * @param version
	 */
	public static void setEcodDatabase(String version) {
		getEcodDatabase(version);
		defaultVersion = version;
	}

	/** Can't instantiate */
	private EcodFactory() {}

}
