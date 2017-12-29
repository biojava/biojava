
/*
 * This file is originally coming from the Dasobert library.
 * Author: Andreas Prlic
 *
 *                  BioJava development code
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
 * Created on Jul 25, 2006
 *
 */

package org.biojava.nbio.structure.align.util;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.net.URL;
import java.net.URLConnection;
import java.util.zip.GZIPInputStream;



/** 
 * A class that takes care about opening URLConnections and sets the proper timeouts
 * @author Andreas Prlic
 * @author Anthony Bradley
 * @since 5.0	
 */
public class URLConnectionTools {

	/** The default connection timeout in ms - 15 seconds*/
	public static final int    DEFAULT_CONNECTION_TIMEOUT = 30000;

	/**
	 * Open HttpURLConnection. Recommended way to open URL connections in Java 1.7 and 1.8.
	 * https://eventuallyconsistent.net/2011/08/02/working-with-urlconnection-and-timeouts/
	 * @param url URL to open
	 * @param timeout timeout in milli seconds
	 * @throws IOException an error in opening the URL
	 */
	public static URLConnection openURLConnection(URL url, int timeout) throws IOException {
		URLConnection huc = url.openConnection();
		huc.setReadTimeout(timeout);
		huc.setConnectTimeout(timeout);
		return huc;
	}


	/** 
	 * Open HttpURLConnection. Recommended way to open
	 * HttpURLConnections, since this take care of setting timeouts
	 * properly for java 1.4 and 1.5
	 * uses the DEFAULT_CONNECTION_TIMEOUT (= 15 seconds)
	 * @param url a URL to open a http connection to
	 * @return HttpURLConnect the opened connection
	 * @throws IOException an error in opening the URL
	 *
	 */
	public static URLConnection openURLConnection(URL url) throws IOException {
		return openURLConnection(url,DEFAULT_CONNECTION_TIMEOUT);
	}

	/** 
	 * Connect to server and return result as an InputStream.
	 * always asks for response to be in GZIP encoded
	 * <p>
	 * The caller is responsible to close the returned InputStream not to cause
	 * resource leaks. 
	 * @param url the URL to connect to
	 * @param timeout the timeout for the connection
	 * @return an {@link InputStream} of response
	 * @throws IOException due to an error opening the URL
	 *
	 */
	public static InputStream getInputStream(URL url, int timeout) throws IOException
	{
		return getInputStream(url,true, timeout);
	}


	/** 
	 * Connect to a URL and return result as an InputStream.
	 * always asks for response to be in GZIP encoded
	 * <p>
	 * The caller is responsible to close the returned InputStream not to cause
	 * resource leaks. 
	 * @param url the input URL to be read
	 * @return an {@link InputStream} of response
	 * @throws IOException due to an error opening the URL
	 */
	public static InputStream getInputStream(URL url) throws IOException
	{
		return getInputStream(url,true, DEFAULT_CONNECTION_TIMEOUT);
	}

	/** 
	 * Open a URL and return an InputStream to it
	 * if acceptGzipEncoding == true, use GZIPEncoding to
	 * compress communication.
	 * <p>
	 * The caller is responsible to close the returned InputStream not to cause
	 * resource leaks.
	 * @param url the input URL to be read
	 * @param acceptGzipEncoding whether to accept Gzip encoding
	 * @param timeout
	 * @return an {@link InputStream} of response
	 * @throws IOException due to an error opening the URL
	 */
	public static InputStream getInputStream(URL url, boolean acceptGzipEncoding, int timeout) throws IOException {
		InputStream inStream = null ;
		URLConnection huc = URLConnectionTools.openURLConnection(url,timeout);

		if ( acceptGzipEncoding) huc.setRequestProperty("Accept-Encoding", "gzip");

		String contentEncoding = huc.getContentEncoding();

		inStream = huc.getInputStream();

		if (contentEncoding != null) {
			if (contentEncoding.contains("gzip")) {
				inStream = new GZIPInputStream(inStream);
			}
		}

		return inStream;

	}

	/** 
	 * Do a POST to a URL and return the response stream for further processing elsewhere.
	 * <p>
	 * The caller is responsible to close the returned InputStream not to cause
	 * resource leaks.  
	 * @param url  the input URL to be read
	 * @param data the post data
	 * @return an {@link InputStream} of response
	 * @throws IOException due to an error opening the URL
	 */
	public static InputStream doPOST(URL url, String data) throws IOException
	{
		return doPOST(url,data,DEFAULT_CONNECTION_TIMEOUT);
	}

	/** 
	 * Do a POST to a URL and return the response stream for further processing elsewhere.
	 * <p>
	 * The caller is responsible to close the returned InputStream not to cause
	 * resource leaks. 
	 * @param url the input URL to be read
	 * @param data the post data
	 * @param timeout
	 * @return an {@link InputStream} of response
	 * @throws IOException due to an error opening the URL
	 */
	public static InputStream doPOST(URL url, String data, int timeout) throws IOException
	{
		URLConnection conn = openURLConnection(url, timeout);
		conn.setDoOutput(true);
		OutputStreamWriter wr = new OutputStreamWriter(conn.getOutputStream());
		wr.write(data);
		wr.flush();
		return conn.getInputStream();
	}


}
