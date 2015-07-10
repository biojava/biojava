package org.biojava.nbio.ws.uniprot;

/**
 * String constants used by this ws
 * @author pbansal
 *
 */
public class UniprotConstants {
	public static final String url = "http://www.uniprot.org";
	public static final String queryUrl = String.format("%s/uniprot/?query=job:%%s&format=%%s", url);
	public static final String queryUrlWithIsoforms = String.format("%s/?query=job:%%s&format=%%s&include=yes", url);
	public static final String uploadListURL = String.format("%s/uploadlists/", url);
} // UniprotConstants
