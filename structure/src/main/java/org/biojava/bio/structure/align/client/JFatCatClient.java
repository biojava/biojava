package org.biojava.bio.structure.align.client;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;


import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLEncoder;

import java.util.Random;
import java.util.SortedSet;
import java.util.TreeSet;


import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.HTTPConnectionTools;
import org.biojava.bio.structure.align.util.ResourceManager;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.bio.structure.align.xml.AFPChainXMLParser;
import org.biojava.bio.structure.align.xml.PdbPairXMLConverter;
import org.biojava.bio.structure.align.xml.RepresentativeXMLConverter;

public class JFatCatClient {

	private static ResourceManager resourceManager = ResourceManager.getResourceManager("jfatcat");

	private static final String serverAPPEND    = "show?name1=%s&name2=%s";
	private static final String sendAPPEND      = "submit?name1=%s&name2=%s&version=%s";
	private static final String multiSendAPPEND = "jobSubmit?username=%s&version=%s";
	private static final String representAPPEND = "representatives?cluster=%s";
	
	
	private static Random generator;
	
	static {
		generator = new Random();
	}
	
	public static AFPChain getAFPChainFromServer(String serverLocation , String name1, String name2, Atom[] ca1, Atom[] ca2) 
	{

	   String serverURL = serverLocation + serverAPPEND;
	   
		try {
			String u = String.format(serverURL,name1,name2);

			URL url = new URL(u);
			System.out.println("requesting alignment from server..."  + url);
			// have a short timeout for this...
			// 5 sec
			InputStream stream = HTTPConnectionTools.getInputStream(url,5000);

			String xml = null;

			if ( stream != null) {

				xml = convertStreamToString(stream);
				//System.out.println("got XML from server: " + xml);
			}
			if (xml != null) {

				//System.out.println("got XML from server: " + xml);

				AFPChain newChain = AFPChainXMLParser.fromXML (xml, ca1, ca2);

				return newChain;

			} else {
				return null;
			} 
		} catch (Exception e){
			System.err.println("error in JFatCatClient: getAFPChainFromServer : " + e.getMessage());
		}
		return null;
	}


	public static String convertStreamToString(InputStream stream){
		BufferedReader reader = new BufferedReader(new InputStreamReader(stream));
		StringBuilder sb = new StringBuilder();

		String line = null;
		try {
			while ((line = reader.readLine()) != null) {
				sb.append(line + "\n");
			}
		} catch (IOException e) {
			//e.printStackTrace();
		} finally {
			try {
				stream.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		return sb.toString();
	}

	public static String sendMultiAFPChainToServer(String serverLocation, String multiXML, String username) throws JobKillException{

	   String multiSendURL = serverLocation + multiSendAPPEND;
	   
		String responseS = "";
		String version = resourceManager.getString("jfatcat.version");

		String u = String.format(multiSendURL,username,version);

		int timeout = getTimeout();

		boolean submitted = false;

		while (! submitted ){
			try { 
				URL url = new URL(u); 
				//System.out.println("posting xml: " + xml);
				InputStream response = HTTPConnectionTools.doPOST(url, multiXML,timeout);
				//System.out.println("got response: " + convertStreamToString(response));
				responseS = convertStreamToString(response);
				submitted = true;
				if (! responseS.contains("OK"))
					System.err.println("server returned " + responseS);

			} catch (Exception e){
				System.err.println("Error in JFatCatClient: while sending results back to server : " + e.getMessage());
				
				try {
					int randomSleep = getRandomSleepTime();
					System.err.println("sleeping " + (randomSleep/1000) + " sec.");
					Thread.sleep(randomSleep);
				} catch (InterruptedException ex){
					ex.printStackTrace();
				}
			}
		} 

		if ( responseS.startsWith("KILL_JOB")){
			throw new JobKillException("Server responded with KILL message.");

		}
		return responseS;
	}

	public static int getRandomSleepTime() {
		
		// we wait between 15 sec and 2 min. 
		
		int minTime = 15000;
		
		int maxTime = 120000 - minTime;
		
		int nextId = generator.nextInt(maxTime);
		return minTime + nextId;
		
	}


	public static final void sendAFPChainToServer(String serverLocation, AFPChain afpChain,Atom[] ca1, Atom[] ca2) 
	{
	   
	   String sendURL = serverLocation + sendAPPEND;

		String version = resourceManager.getString("jfatcat.version");
		try {

			// just to make sure that similarity has been calculated!
			afpChain.getSimilarity();

			String xml = AFPChainXMLConverter.toXML(afpChain, ca1, ca2);

			String u = String.format(sendURL,afpChain.getName1() , afpChain.getName2(),version);

			URL url = new URL(u); 
			//System.out.println("posting xml: " + xml);
			InputStream response = HTTPConnectionTools.doPOST(url, xml);
			FarmJobRunnable.log("got response: " + convertStreamToString(response));
			if ( xml.startsWith("KILL_JOB")){
				throw new JobKillException("Server responded with KILL message.");
			}

		} catch (Exception e){
			System.err.println("error in JFatCatClient: sendAFPChainToServer : " + e.getMessage());
		}

	}

	public static final int getTimeout(){
		String timeoutS = resourceManager.getString("connection.timeout");
		int timeout = 60000;

		try {
			timeout = Integer.parseInt(timeoutS);
		} catch (NumberFormatException ex ){
			ex.printStackTrace();
		}
		return timeout;
	}

	public static final SortedSet<PdbPair> getPdbPairs(String url,int nrPair, String username) throws MalformedURLException, IOException, JobKillException {
		StringBuffer u = new StringBuffer();
		u.append(url);
		u.append("getPairs?");
		u.append("nrPairs=");
		u.append(nrPair);
		u.append("&username=");
		u.append(URLEncoder.encode(username,"UTF-8"));


		String urlS= u.toString();
		int timeout = getTimeout();

		SortedSet<PdbPair> pairs = new TreeSet<PdbPair>();
		FarmJobRunnable.log("requesting " + urlS);
		URL serverUrl = new URL(urlS);
		// we are very tolerant with requesting a set of pairs, since we probably depend on getting it to get work started...
		// 1 min...
		InputStream stream = HTTPConnectionTools.getInputStream(serverUrl,timeout);
		String xml = null;

		if ( stream != null) {

			xml = convertStreamToString(stream);
			//System.out.println("got XML from server: " + xml);
			if (xml != null) {
				if ( xml.startsWith("KILL_JOB")){
					// we got the KILL signal from the server...
					throw new JobKillException("Server responded with KILL message.");
				}
				pairs = PdbPairXMLConverter.convertXMLtoPairs(xml);
			} 
		}

		return pairs;
	}


	public static final SortedSet<String> getRepresentatives(String serverLocation, int cutoff){
		SortedSet<String> representatives = new TreeSet<String>();
		
		String representURL = serverLocation + representAPPEND;
		
		if ( cutoff < 20)
			cutoff = 40;
		int timeout = getTimeout();
		String u = String.format(representURL,cutoff);
		try {
			URL url = new URL(u);
			
			InputStream stream = HTTPConnectionTools.getInputStream(url,timeout);

			String xml = null;

			if ( stream != null) {

				xml = convertStreamToString(stream);
				//System.out.println("got XML from server: " + xml);
			}
			if (xml != null) {
				representatives = RepresentativeXMLConverter.fromXML(xml);
			}
		} catch (Exception e){
			e.printStackTrace();
		}
		return representatives;
	}



}
