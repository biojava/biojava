/*
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
package org.biojava.dasobert.util;

import java.io.IOException;
import java.io.InputStream;
import java.lang.reflect.Method;
import java.net.ConnectException;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.zip.GZIPInputStream;

import org.biojava.bio.program.das.dasalignment.DASException;
import org.biojava.dasobert.das.DASStatus;


/** a class that takes care about opening HttpURLConnections and sets the proper timeouts
 * 
 * @author Andreas Prlic
 * @since 9:58:25 AM
 * @version %I% %G%
 */
public class HttpConnectionTools {



	static int    DEFAULT_CONNECTION_TIMEOUT = 15000; // timeout for http connection = 15. sec


	public HttpConnectionTools() {
		super();

	}

	/**open HttpURLConnection. Recommended way to open
	 * HttpURLConnections, since this take care of setting timeouts
	 * properly for java 1.4 and 1.5
	 * 
	 * @param url URL to oopen
	 * @param timeout timeout in milli seconds
	 * @return a HttpURLConnection
	 * @throws IOException 
	 * @throws ConnectException 
	 * 
	 *
	 */
	public static HttpURLConnection openHttpURLConnection(URL url, int timeout)
	throws IOException, ConnectException{

		HttpURLConnection huc = null;

		huc = (HttpURLConnection) url.openConnection();
		huc.addRequestProperty("User-Agent", DasobertDefaults.USERAGENT);

		// this sets the timeouts for Java 1.4
		System.setProperty("sun.net.client.defaultConnectTimeout", ""+timeout);
		System.setProperty("sun.net.client.defaultReadTimeout", ""+timeout);

		// for Java 1.5 we need to do this:
		// use reflection to determine if get and set timeout methods for urlconnection are available
		// seems java 1.5 does not watch the System properties any longer...
		// and java 1.4 did not provide the new classes...

		try {
			// try to use reflection to set timeout property
			Class urlconnectionClass = Class.forName("java.net.HttpURLConnection");

			Method setconnecttimeout = urlconnectionClass.getMethod (
					"setConnectTimeout", new Class [] {int.class}        
			); 
			setconnecttimeout.invoke(huc,new Object[] {new Integer(timeout)});

			Method setreadtimeout = urlconnectionClass.getMethod (
					"setReadTimeout", new Class[] {int.class}
			);
			setreadtimeout.invoke(huc,new Object[] {new Integer(timeout)});
			//System.out.println("successfully set java 1.5 timeout");
		} catch (Exception e) {
			e.printStackTrace();
			// most likely it was a NoSuchMEthodException and we are running java 1.4.
		}
		return huc;
	}


	/** open HttpURLConnection. Recommended way to open
	 * HttpURLConnections, since this take care of setting timeouts
	 * properly for java 1.4 and 1.5
	 * uses the DEFAULT_CONNECTION_TIMEOUT (= 15 seconds)
	 * 
	 * @param url a URL to open a http connection to
	 * @return HttpURLConnect the opened connection
	 * @throws IOException
	 * @throws ConnectException
	 * 
	 * */
	public static HttpURLConnection openHttpURLConnection(URL url) 
	throws IOException, ConnectException {

		return openHttpURLConnection(url,DEFAULT_CONNECTION_TIMEOUT);

	}

	/** connect to DAS server and return result as an InputStream.
	 * always asks for response to be in GZIP encoded
	 * 
	 * @param url the URL to connect to
	 * @return an InputStream
	 * @throws IOException 
	 * @throws DASException if DAS server returns error response code
	 *
	 */    
	public static InputStream getInputStream(URL url) 
	throws IOException, DASException
	{
		return getInputStream(url,true);
	}

	/** open a URL and return an InputStream to it
	 *  if acceptGzipEncoding == true, use GZIPEncoding to
	 *  compress communication
	 * 
	 * @param url
	 * @param acceptGzipEncoding
	 * @return an InputStream to the URL
	 * @throws IOException
	 * @throws DASException if DAS server returns error response code 
	 */
	public static InputStream getInputStream(URL url, boolean acceptGzipEncoding)
	throws IOException, DASException {
		InputStream inStream = null ;

		//System.out.println("opening connection to "+url);
		HttpURLConnection huc = HttpConnectionTools.openHttpURLConnection(url);  

		if ( acceptGzipEncoding) {
			// should make communication faster
			huc.setRequestProperty("Accept-Encoding", "gzip");
		}


		// check the HTTP response code
		int httpCode = huc.getResponseCode();

		// and the DAS response code which is set in the
		// http header X-DAS-Status       
		String dasCodeS = huc.getHeaderField("X-DAS-Status");
		if ( dasCodeS != null) {
			int dasCode = -1;
			try {
				dasCode = Integer.parseInt(dasCodeS);
			} catch (NumberFormatException e){
				// ignroe it
			}
			String desc = DASStatus.getErrorDescription(dasCode);

			if ( ! desc.equals(DASStatus.UNKNOWN))
				if (( httpCode == HttpURLConnection.HTTP_OK) && 
						( dasCode != DASStatus.STATUS_OK)){

					throw new DASException(desc);

				}
		}

		String contentEncoding = huc.getContentEncoding();

		inStream = huc.getInputStream();

		if (contentEncoding != null) {
			if (contentEncoding.indexOf("gzip") != -1) {
				// we have gzip encoding
				inStream = new GZIPInputStream(inStream);               
			}
		}

		return inStream;

	}


	/** request an InputStream and provide username and password for logging in.
	 *  (using BASE64 Encoding)
	 * 
	 * @param url
	 * @param name
	 * @param password
	 * @return an InputStream
	 * @throws IOException
	 * @throws ConnectException
	 */
	public static InputStream getInputStream(URL url, String name, String password)
	throws IOException, ConnectException
	{
		InputStream inStream = null;


		HttpURLConnection huc = openHttpURLConnection(url);

		huc.setRequestProperty(
				"Authorization", 
				"Basic " + 
				encode(name + ":" + password)
		);

		inStream = huc.getInputStream();        

		return inStream;

	}



	public static String encode (String source) {
		
		return(BASE64Encoder.encodeString(source));
	}
}
