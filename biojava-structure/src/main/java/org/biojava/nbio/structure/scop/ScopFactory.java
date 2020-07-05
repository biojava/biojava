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
package org.biojava.nbio.structure.scop;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.HashMap;
import java.util.Map;



/**
 * Controls the global ScopDatabase being used.
 *
 * <p>
 * Example: Fetch the structure corresponding to an old version of scop
 *
 * <pre>
 * ScopInstallation scop = new ScopInstallation();
 * scop.setScopVersion("1.69");
 * ScopFactory.setScopDatabase(scop);
 * AtomCache cache = new AtomCache();
 * cache.setFetchFileEvenIfObsolete(true); //fetch older PDBs
 * cache.setStrictSCOP(false); // correct simple errors in domain names
 * Structure s = cache.getStructure("d3hbia_");
 * </pre>
 * @author sbliven
 *
 */
public class ScopFactory {

	private static final Logger logger = LoggerFactory.getLogger(ScopFactory.class);

	// berkeley 2
	public static final String VERSION_2_0_7 = "2.07";
	public static final String VERSION_2_0_6 = "2.06";
	public static final String VERSION_2_0_5 = "2.05";
	public static final String VERSION_2_0_4 = "2.04";
	public static final String VERSION_2_0_3 = "2.03";
	public static final String VERSION_2_0_2 = "2.02";
	public static final String VERSION_2_0_1 = "2.01";
	// berkeley 1 (aliases of above)
	public static final String VERSION_1_75C = VERSION_2_0_3;
	public static final String VERSION_1_75B = VERSION_2_0_2;
	public static final String VERSION_1_75A = VERSION_2_0_1;
	// original SCOP
	// latest SCOP release from SCOP website = 1.75;
	public static final String VERSION_1_75 = "1.75";
	// outdated stable versions
	public static final String VERSION_1_73 = "1.73";
	public static final String VERSION_1_71 = "1.71";
	public static final String VERSION_1_69 = "1.69";
	public static final String VERSION_1_67 = "1.67";
	public static final String VERSION_1_65 = "1.65";
	public static final String VERSION_1_63 = "1.63";
	public static final String VERSION_1_61 = "1.61";
	public static final String VERSION_1_59 = "1.59";
	public static final String VERSION_1_57 = "1.57";
	public static final String VERSION_1_55 = "1.55";

	// The most recent version as of compilation time
	public static final String LATEST_VERSION = VERSION_2_0_7;

	// Hold one instance for each version
	private static Map<String,ScopDatabase> versionedScopDBs = new HashMap<>();
	private static String defaultVersion = LATEST_VERSION;

	/**
	 * Get the current default instance for the default version
	 * @return
	 */
	public static ScopDatabase getSCOP(){
		return getSCOP(defaultVersion);
	}

	/**
	 * Gets an instance of the specified scop version.
	 *
	 * <p>
	 * The particular implementation returned is guaranteed to
	 * implement {@link LocalScopDatabase} (generally a {@link BerkeleyScopInstallation}).
	 *
	 * @param version A version number, such as {@link #VERSION_1_75A}
	 * @return an
	 */
	public static ScopDatabase getSCOP(String version){
		if( version == null ) {
			version = defaultVersion;
		}

		ScopDatabase scop = versionedScopDBs.get(version);
		if (scop == null) {
			// Use a local installation
			logger.info("Creating new {}, version {}", BerkeleyScopInstallation.class.getSimpleName(), version);
			BerkeleyScopInstallation berkeley = new BerkeleyScopInstallation();
			berkeley.setScopVersion(version);
			versionedScopDBs.put(version, berkeley);
			return berkeley;
		}
		return scop;
	}

	/**
	 * Set the default scop version
	 * @param version A version number, such as {@link #VERSION_1_75A}
	 */
	public static void setScopDatabase(String version) {
		getSCOP(version);
		defaultVersion = version;
	}

	/**
	 * Set the default scop version and instance
	 * @param scop
	 */
	public static void setScopDatabase(ScopDatabase scop){
		logger.debug("ScopFactory: Setting ScopDatabase to type: {}", scop.getClass().getName());
		defaultVersion = scop.getScopVersion();
		versionedScopDBs.put(defaultVersion,scop);
	}
}
