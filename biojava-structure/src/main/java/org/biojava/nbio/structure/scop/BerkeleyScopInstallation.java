/**
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
 * Created on Oct 12, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.nbio.structure.scop;

import java.util.HashMap;
import java.util.Map;

/**
 * SCOPe:
 *
 * The Structural Classification of Proteins (extended) at Berkeley Lab and UC Berkeley
 * (<a href="http://scop.berkeley.edu/">http://scop.berkeley.edu/</a>).
 *
 * <p>This provides updates to the MRC SCOP hierarchical classification.
 */
public class BerkeleyScopInstallation extends ScopInstallation {


	String defaultBerkeleyDownloadURL = "http://scop.berkeley.edu/downloads/parse/";
	String defaultBerkeleyScopVersion=ScopFactory.LATEST_VERSION;

	/**
	 * A map from SCOP version names which the Berkeley server offers as a
	 * download to an array of equivalent deprecated SCOP version names.
	 */
	public static final Map<String,String[]> EQUIVALENT_VERSIONS = new HashMap<String,String[]>();

	static {
		EQUIVALENT_VERSIONS.put("2.01", new String[] {"1.75A"});
		EQUIVALENT_VERSIONS.put("2.02", new String[] {"1.75B"});
		EQUIVALENT_VERSIONS.put("2.03", new String[] {"1.75C"});
	}


	public BerkeleyScopInstallation() {
		super();
		setScopVersion(defaultBerkeleyScopVersion);
		addMirror(new BerkeleyScopMirror(defaultBerkeleyDownloadURL));
		addMirror(new ScopMirror());
	}

	private static class BerkeleyScopMirror extends ScopMirror {
		private String rootURL;
		public BerkeleyScopMirror(String url) {
			super(url);
			rootURL = url;
		}

		@Override
		public String getClaURL(String scopVersion) {
			return rootURL+getFilename("cla",scopVersion);
		}

		@Override
		public String getDesURL(String scopVersion) {
			return rootURL+getFilename("des",scopVersion);
		}

		@Override
		public String getHieURL(String scopVersion) {
			return rootURL+getFilename("hie",scopVersion);
		}

		@Override
		public String getComURL(String scopVersion) {
			return rootURL+getFilename("com",scopVersion);
		}

		private String getFilename(String fileType, String version) {
			// Convert to canonical version number
			for (Map.Entry<String, String[]> entry : EQUIVALENT_VERSIONS.entrySet()) {
				for (String vr : entry.getValue()) {
					if (version.equals(vr)) {
						version = entry.getKey();
						break;
					}
				}
			}
			String[] parts = version.split("\\.");
			// they changed the filename schemes!
			if (Integer.parseInt(parts[0]) == 1) {
				return "dir." + fileType + ".scop." + version + ".txt";
			} else {
				return "dir." + fileType + ".scope." + version + "-stable.txt";
			}
		}

		@Override
		public String toString() {
			return "BerkeleyScopMirror[ \"" + rootURL + " ]";
		}
	}

}
