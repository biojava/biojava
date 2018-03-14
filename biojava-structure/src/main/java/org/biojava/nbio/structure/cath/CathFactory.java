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
package org.biojava.nbio.structure.cath;

import org.biojava.nbio.structure.scop.ScopFactory;

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

	public static final String VERSION_3_5_0 = "3_5_0";
	public static final String VERSION_4_0_0 = "4_0_0";
	public static final String VERSION_4_1_0 = "4_1_0";
	public static final String VERSION_4_2_0 = "4_2_0";
	public static final String LATEST_VERSION = VERSION_4_2_0;

	public static final String DEFAULT_VERSION = LATEST_VERSION;

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
