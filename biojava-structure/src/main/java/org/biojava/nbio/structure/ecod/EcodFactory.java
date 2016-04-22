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
import java.lang.ref.SoftReference;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;

import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.cath.CathFactory;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Controls global {@link EcodDatabase EcodDatabases} being used.
 * Implements a multiton pattern through {@link #getEcodDatabase(String)},
 * and a singleton pattern through {@link #getEcodDatabase()}.
 * @author Spencer Bliven
 * @see ScopFactory
 * @see CathFactory
 * @see EcodInstallation
 */
public class EcodFactory {

	private static final Logger logger = LoggerFactory.getLogger(EcodFactory.class);
	
	public static final String DEFAULT_VERSION = EcodInstallation.DEFAULT_VERSION;

	private static Map<String, SoftReference<EcodDatabase>> versionedEcodDBs =
			Collections.synchronizedMap(new HashMap<String, SoftReference<EcodDatabase>>());
	private static String defaultVersion = EcodInstallation.DEFAULT_VERSION;

	/**
	 * Returns the (singleton) database for the current default version
	 */
	public static EcodDatabase getEcodDatabase() {
		return getEcodDatabase(defaultVersion);
	}

	public static EcodDatabase getEcodDatabase(String version) {
		if( version == null )
			version = defaultVersion;

		logger.trace("Waiting for EcodFactory lock to get version "+version);
		synchronized(versionedEcodDBs) {
			logger.trace("Got EcodFactory lock to get version "+version);

			releaseReferences();

			SoftReference<EcodDatabase> ecodRef = versionedEcodDBs.get(version.toLowerCase());
			EcodDatabase ecod = null;
			if(ecodRef != null) {
				ecod = ecodRef.get();
			}
			if( ecod == null ) {
				logger.debug("Creating new {}, version {}",EcodInstallation.class.getSimpleName(),version);
				String cacheDir = new UserConfiguration().getCacheFilePath();
				ecod = new EcodInstallation(cacheDir, version);
				versionedEcodDBs.put(version.toLowerCase(), new SoftReference<EcodDatabase>(ecod));

				// If the parsed version differed from that requested, add that too
				// Note that getVersion() may trigger file parsing
				try {
					if( ! versionedEcodDBs.containsKey(ecod.getVersion().toLowerCase()) ) {
						versionedEcodDBs.put(ecod.getVersion().toLowerCase(),new SoftReference<EcodDatabase>(ecod));
					}
				} catch (IOException e) {
					// For parsing errors, just use the requested version
				}
			}
			logger.trace("Releasing EcodFactory lock after getting version "+version);

			return ecod;
		}
	}

	/**
	 * removes SoftReferences which have already been garbage collected
	 */
	private static void releaseReferences() {
		synchronized(versionedEcodDBs) {
			Iterator<Entry<String, SoftReference<EcodDatabase>>> it = versionedEcodDBs.entrySet().iterator();
			while(it.hasNext()) {
				Entry<String, SoftReference<EcodDatabase>> entry = it.next();
				SoftReference<EcodDatabase> ref = entry.getValue();
				if(ref.get() == null) {
					logger.debug("Removed version {} from EcodFactory to save memory.",entry.getKey());
					it.remove();
				}
			}
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
