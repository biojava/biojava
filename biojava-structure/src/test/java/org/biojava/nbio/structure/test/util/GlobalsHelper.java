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

import java.util.Deque;
import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;

import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.ChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava.nbio.structure.scop.ScopDatabase;
import org.biojava.nbio.structure.scop.ScopFactory;

/**
 * Helper class to manage all the global state changes in BioJava.
 * For instance, this should be used in tests before modifying PDB_PATH.
 *
 * Used by tests during setup and teardown to ensure a clean environment
 *
 * This class is a singleton.
 * @author Spencer Bliven
 *
 */
public final class GlobalsHelper {

	private static class PathInfo {
		public final String pdbPath;
		public final String pdbCachePath;
		public final AtomCache atomCache;
		public final ChemCompProvider chemCompProvider;
		public final String downloadChemCompProviderPath;
		public final ScopDatabase scop;

		public PathInfo() {
			pdbPath = System.getProperty(UserConfiguration.PDB_DIR, null);
			pdbCachePath = System.getProperty(UserConfiguration.PDB_CACHE_DIR, null);
			atomCache = StructureIO.getAtomCache();
			chemCompProvider = ChemCompGroupFactory.getChemCompProvider();
			downloadChemCompProviderPath = DownloadChemCompProvider.getPath().getPath();
			scop = ScopFactory.getSCOP();
		}
	}

	// Saves defaults as stack
	private static Deque<PathInfo> stack = new LinkedList<>();
	static {
		// Save default state
		pushState();
	}

	/**
	 * GlobalsHelper should not be instantiated.
	 */
	private GlobalsHelper() {}

	/**
	 * Save current global state to the stack
	 */
	public static void pushState() {
		PathInfo paths = new PathInfo();
		stack.addFirst(paths);
	}

	/**
	 * Sets a new PDB_PATH and PDB_CACHE_PATH consistently.
	 *
	 * Previous values can be restored with {@link #restoreState()}.
	 * @param path
	 */
	public static void setPdbPath(String path, String cachePath) {
		pushState();
		if(path == null || cachePath == null) {
			UserConfiguration config = new UserConfiguration();
			if(path == null) {
				path = config.getPdbFilePath();
			}
			if(cachePath == null) {
				cachePath = config.getCacheFilePath();
			}
		}
		System.setProperty(UserConfiguration.PDB_DIR, path);
		System.setProperty(UserConfiguration.PDB_CACHE_DIR, path);

		AtomCache cache = new AtomCache(path);
		StructureIO.setAtomCache(cache);

		// Note side effect setting the path for all DownloadChemCompProvider due to static state
		ChemCompProvider provider = new DownloadChemCompProvider(path);
		ChemCompGroupFactory.setChemCompProvider(provider);
	}

	/**
	 * Restore global state to the previous settings
	 * @throws NoSuchElementException if there is no prior state to restore
	 */
	public static void restoreState() {
		PathInfo paths = stack.removeFirst();

		if(paths.pdbPath == null) {
			System.clearProperty(UserConfiguration.PDB_DIR);
		} else {
			System.setProperty(UserConfiguration.PDB_DIR, paths.pdbPath);
		}
		if(paths.pdbCachePath == null) {
			System.clearProperty(UserConfiguration.PDB_CACHE_DIR);
		} else {
			System.setProperty(UserConfiguration.PDB_CACHE_DIR, paths.pdbCachePath);
		}

		StructureIO.setAtomCache(paths.atomCache);

		// Use side effect setting the path for all DownloadChemCompProvider due to static state
		new DownloadChemCompProvider(paths.downloadChemCompProviderPath);

		ChemCompGroupFactory.setChemCompProvider(paths.chemCompProvider);

		ScopFactory.setScopDatabase(paths.scop);
	}


}
