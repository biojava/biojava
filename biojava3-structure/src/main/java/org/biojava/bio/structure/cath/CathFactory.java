package org.biojava.bio.structure.cath;

import java.util.HashMap;
import java.util.Map;

/**
 * Controls global {@link CathDatabase CathDatabases} being used.
 * Implements a multiton pattern through {@link #getCathDatabase(String)},
 * and a singleton pattern through {@link #getCathDatabase()}.
 * @author dmyersturnbull
 * @see ScopFactory
 * @see CathInstallation
 */
public class CathFactory {

	public static final String VERSION_3_5_0 = "3.5.0";
	public static final String LATEST_VERSION = VERSION_3_5_0;
	
	public static String DEFAULT_VERSION = LATEST_VERSION;
	
	private static CathDatabase cath;
	
	private static Map<String, CathDatabase> versions = new HashMap<String, CathDatabase>();
	
	/**
	 * Sets the default (singleton) CathDatabase.
	 */
	public static void setCath(CathDatabase cath) {
		CathFactory.cath = cath;
	}

	/**
	 * Returns the default (singleton) CathDatabase.
	 * If the database is null, this will recreate it (lazy initialization).
	 */
	public static CathDatabase getCathDatabase() {
		if (cath == null) {
			cath = new CathInstallation();
		}
		return cath;
	}
	
	private CathFactory() {
		
	}

	/**
	 * Returns a CATH database of the specified version.
	 * @param version For example, "3.5.0"
	 */
	public static CathDatabase getCathDatabase(String version) {
		if (version == null) version = DEFAULT_VERSION;
		CathDatabase cath = versions.get(version);
		if (cath == null) {
			CathInstallation newCath = new CathInstallation();
			newCath.setCathVersion(version);
			cath = newCath;
		}
		return cath;
	}
	
}
