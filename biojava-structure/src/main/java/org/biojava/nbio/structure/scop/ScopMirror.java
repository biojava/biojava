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

import org.biojava.nbio.core.util.FileDownloadUtils;

import java.net.URI;
import java.net.URISyntaxException;

/**
 * Helper class to store paths to the four SCOP files
 *
 * The string "%s" is replaced with the version number.
 * @author Spencer Bliven
 *
 */
public class ScopMirror {
	private String rootURL;
	private final String claURL;
	private final String desURL;
	private final String hieURL;
	private final String comURL;

	/** Specify all keys individually */
	public ScopMirror(String claURL, String desURL,
			String hieURL, String comURL) {
		this.rootURL = null;
		this.claURL = claURL;
		this.desURL = desURL;
		this.hieURL = hieURL;
		this.comURL = comURL;
	}
	/** Specify a common root URL which is concatenated with individual filenames */
	public ScopMirror(String url, String claURL, String desURL,
			String hieURL, String comURL) {
		this(url+claURL, url+desURL, url+hieURL, url+comURL);
		this.rootURL = url;
	}
	public ScopMirror(String url) {
		this(url,
				"dir.cla.scop.txt_%s","dir.des.scop.txt_%s",
				"dir.hie.scop.txt_%s","dir.com.scop.txt_%s");
	}
	/** Use default MRC location */
	public ScopMirror() {
		this(ScopInstallation.SCOP_DOWNLOAD);
	}

	/**
	 * Get the URL for the root download directory, or null if none is set.
	 * @return
	 */
	public String getRootURL() {
		return rootURL;
	}
	public String getClaURL(String scopVersion) {
		return String.format(claURL,scopVersion);
	}
	public String getDesURL(String scopVersion) {
		return String.format(desURL,scopVersion);
	}
	public String getHieURL(String scopVersion) {
		return String.format(hieURL,scopVersion);
	}
	public String getComURL(String scopVersion) {
		return String.format(comURL,scopVersion);
	}
	public boolean isReachable() {
		final int timeout = 800;
		if(rootURL != null) {
			return FileDownloadUtils.ping(getRootURL(),timeout);
		} else {
			try {
				URI cla = new URI(getClaURL("VERSION"));
				String host = cla.getHost();
				return FileDownloadUtils.ping(host,timeout);
			} catch (URISyntaxException e) {
				return false;
			}
		}
	}
}
