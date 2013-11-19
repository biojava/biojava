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
package org.biojava.bio.structure.scop;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

/**
 * The Structural Classification of Proteins at Berkeley Lab and UC Berkeley (<a href="http://scop.berkeley.edu/">http://scop.berkeley.edu/</a>).
 */
public class BerkeleyScopInstallation extends ScopInstallation {


	String defaultBerkeleyDownloadURL = "http://scop.berkeley.edu/downloads/parse/";
	String defaultBerkeleyScopVersion="2.0.3";

	/**
	 * A map from SCOP version names which the Berkeley server offers as a download to an array of equivalent deprecated SCOP version names.
	 */
	public static final Map<String,String[]> EQUIVALENT_VERSIONS = new HashMap<String,String[]>();

	static {
		EQUIVALENT_VERSIONS.put("2.0.1", new String[] {"1.75A"});
		EQUIVALENT_VERSIONS.put("2.0.2", new String[] {"1.75B"});
		EQUIVALENT_VERSIONS.put("2.0.3", new String[] {"1.75C"});
	}


	public BerkeleyScopInstallation() {
		super();
		setScopVersion(defaultBerkeleyScopVersion);
		setScopDownloadURL(defaultBerkeleyDownloadURL);
	}

	private String getFilename(String fileType) {
		String version = scopVersion;
		for (Map.Entry<String, String[]> entry : EQUIVALENT_VERSIONS.entrySet()) {
			for (String vr : entry.getValue()) {
				if (scopVersion.equals(vr)) {
					version = entry.getKey();
					break;
				}
			}
		}
		String[] parts = version.split("\\.");
		if (Integer.parseInt(parts[0]) == 1) {
			return "dir." + fileType + ".scop." + version + ".txt";
		} else {
			return "dir." + fileType + ".scope." + version + "-stable.txt";
		}
	}


	protected void downloadClaFile() throws IOException {
		String filename = getFilename("cla");
		URL url = new URL(scopDownloadURL + filename);
		downloadFileFromRemote(url, new File(filename));
	}

	protected void downloadDesFile() throws IOException {
		String filename = getFilename("des");
		URL url = new URL(scopDownloadURL + filename);
		downloadFileFromRemote(url, new File(filename));
	}

	protected void downloadHieFile() throws IOException {
		String filename = getFilename("hie");
		URL url = new URL(scopDownloadURL + filename);
		downloadFileFromRemote(url, new File(filename));
	}
	
	protected void downloadComFile() throws IOException {
		String filename = getFilename("com");
		URL url = new URL(scopDownloadURL + filename);
		downloadFileFromRemote(url, new File(filename));
	}

}
