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
import org.biojava.bio.structure.align.fatcat.FatCatRigid;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.HTTPConnectionTools;
import org.biojava.bio.structure.align.util.ResourceManager;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.bio.structure.align.xml.AFPChainXMLParser;
import org.biojava.bio.structure.align.xml.HasResultXMLConverter;
import org.biojava.bio.structure.align.xml.PdbPairXMLConverter;
import org.biojava.bio.structure.align.xml.PdbPairsMessage;
import org.biojava.bio.structure.align.xml.PositionInQueueXMLConverter;
import org.biojava.bio.structure.align.xml.RepresentativeXMLConverter;

public class JFatCatClient {

	private static ResourceManager resourceManager = ResourceManager.getResourceManager("jfatcat");

	private static final String serverAPPEND    = "show?name1=%s&name2=%s";
	
	private static final String sendAPPEND      = "submit?name1=%s&name2=%s&version=%s";
	
	private static final String multiSendAPPEND = "jobSubmit?username=%s&version=%s";
	
	private static final String representAPPEND = "representatives?cluster=%s";
	
	private static final String serverHasResult = "hasResult?method=%s&name1=%s&name2=%s";

	private static final int DEFAULT_TIMEOUT = 5000;

	private static final String serverPositionInQueue =  "queuePosition?method=%s&name1=%s&name2=%s";

	private static Random generator;

	private static String newline = System.getProperty("line.separator");

	private static String KILL_JOB = "KILL_JOB";
	
	private static String COME_BACK_LATER = "COME_BACK_LATER";

	static {
		
		generator = new Random();
	
	}

	static final boolean debug = false;

	public static void main(String[] args){
		//System.out.println(hasPrecalculatedResult("http://emmy.rcsb.org/jfatcatserver/align/", "jCE Circular Permutation", "1CDG.A", "1TIM.A"));
		AtomCache cache = new AtomCache();
		String name1= "2W72.A";
		String name2= "1D2Z.D";
		
		try {
		
			Atom[] ca1 = cache.getAtoms(name1);
			Atom[] ca2 = cache.getAtoms(name2);
			
			int timeout = 10000;

			String testServer = "http://emmy.rcsb.org/jfatcatserver/align/";
			
			System.out.println(getAFPChainFromServer(testServer, FatCatRigid.algorithmName, name1, name2, ca1, ca2, timeout));

			PdbPairsMessage msg = getPdbPairs(testServer, 1, "test");
			
			System.out.println(msg);

		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public static boolean hasPrecalculatedResult(String serverLocation, String method, String name1, String name2 ){
		return hasPrecalculatedResult(serverLocation, method, name1, name2, DEFAULT_TIMEOUT );
	}

	public static boolean hasPrecalculatedResult(String serverLocation, String method, String name1, String name2, int timeout){

		String serverURL = serverLocation + serverHasResult;


		boolean hasResults = false;
		try {
			String u = String.format(serverURL,URLEncoder.encode(method,"UTF-8"),name1,name2) ;
			URL url = new URL(u);
			//System.out.println("has result ? ..."  + url);

			InputStream stream = HTTPConnectionTools.getInputStream(url,timeout);

			String xml = null;

			if ( stream != null) {

				xml = convertStreamToString(stream);
				if ( debug )
					System.out.println(" has PrecalcResults got XML from server: " + xml);
				HasResultXMLConverter conv = new HasResultXMLConverter();
				hasResults = conv.fromXML(xml);
			}

		} catch (Exception e){
			System.err.println("error in JFatCatClient: getAFPChainFromServer : " + e.getMessage());
		}
		return hasResults;
	}


	public int getPositionInQueue(String serverLocation, String method, String name1, String name2){
		return getPositionInQueue(serverLocation, method, name1, name2, DEFAULT_TIMEOUT);
	}

	public int getPositionInQueue(String serverLocation, String method, String name1, String name2, int timeout){
		String serverURL = serverLocation + serverPositionInQueue;


		int position = Integer.MIN_VALUE;
		try {
			String u = String.format(serverURL,URLEncoder.encode(method,"UTF-8"),name1,name2) ;
			URL url = new URL(u);

			InputStream stream = HTTPConnectionTools.getInputStream(url,timeout);

			String xml = null;

			if ( stream != null) {

				xml = convertStreamToString(stream);
				//System.out.println("got XML from server: " + xml);
				PositionInQueueXMLConverter conv = new PositionInQueueXMLConverter();
				position = conv.fromXML(xml);
			}

		} catch (Exception e){
			System.err.println("error in JFatCatClient: getAFPChainFromServer : " + e.getMessage());
		}
		return position;

	}
	public static AFPChain getAFPChainFromServer(String serverLocation ,  String name1, String name2, Atom[] ca1, Atom[] ca2) {
		String method = FatCatRigid.algorithmName;
		return getAFPChainFromServer(serverLocation, method, name1, name2, ca1, ca2,DEFAULT_TIMEOUT);
	}

	public static AFPChain getAFPChainFromServer(String serverLocation , String method, String name1, String name2, Atom[] ca1, Atom[] ca2, int timeout) 
	{

		String serverURL = serverLocation + serverAPPEND;

		try {
			String u = String.format(serverURL,name1,name2) ;

			if ( method != null)
				u+= "&method=" + URLEncoder.encode(method,"UTF-8");

			URL url = new URL(u);
			System.out.println("requesting alignment from server..."  + url);
			// have a short timeout for this...
			// 5 sec
			InputStream stream = HTTPConnectionTools.getInputStream(url,timeout);

			String xml = null;

			if ( stream != null) {

				xml = convertStreamToString(stream);
				//System.out.println("got XML from server: " + xml);
			}
			if (xml != null) {

				//System.out.println("got XML from server: " + xml);

				AFPChain newChain = AFPChainXMLParser.fromXML (xml, name1, name2, ca1, ca2);

				return newChain;

			} else {
				return null;
			} 
		} catch (Exception e){
			e.printStackTrace();
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
                sb.append(line).append(newline);
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
		String version = resourceManager.getString("jfatcat.version");
		return sendMultiAFPChainToServer(serverLocation, multiXML, username, version);
	}
	
	public static String sendMultiAFPChainToServer(String serverLocation, String multiXML, String username, String version) throws JobKillException{
		String multiSendURL = serverLocation + multiSendAPPEND;

		String responseS = "";
		
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

				// server is busy... wait a bit and try again
				if ( responseS.startsWith(COME_BACK_LATER)){
					submitted = false;
				}

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

		if ( responseS.startsWith(KILL_JOB)){
			throw new JobKillException("Server responded with KILL message.");

		}


		return responseS;
	}

	public static int getRandomSleepTime() {

		// now wait between 7 and 13 min. 

		int minTime = 560000;

		int maxTime = 7800000 - minTime;

		int nextId = generator.nextInt(maxTime);
		return minTime + nextId;

	}


	public static final void sendAFPChainToServer(String serverLocation, AFPChain afpChain,Atom[] ca1, Atom[] ca2) 
	{

		String sendURL = serverLocation + sendAPPEND;

		String version = resourceManager.getString("jfatcat.version");

		int timeout = getTimeout();

		try {

			// just to make sure that similarity has been calculated!
			afpChain.getSimilarity();

			String xml = AFPChainXMLConverter.toXML(afpChain, ca1, ca2);

			String u = String.format(sendURL,afpChain.getName1() , afpChain.getName2(),version);

			URL url = new URL(u); 

			//System.out.println("posting xml: " + xml);
			InputStream response = HTTPConnectionTools.doPOST(url, xml,timeout);

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


	public static final PdbPairsMessage getPdbPairs(String url,int nrPair, String username) throws MalformedURLException, IOException, JobKillException {
		StringBuffer u = new StringBuffer();
		u.append(url);
		u.append("getPairs?");
		u.append("nrPairs=");
		u.append(nrPair);
		u.append("&username=");
		u.append(URLEncoder.encode(username,"UTF-8"));


		String urlS= u.toString();
		int timeout = getTimeout();

		//SortedSet<PdbPair> pairs = new TreeSet<PdbPair>();
		PdbPairsMessage msg = null;
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
				msg = PdbPairXMLConverter.convertXMLtoPairs(xml);

			} 
		}

		return msg;
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
