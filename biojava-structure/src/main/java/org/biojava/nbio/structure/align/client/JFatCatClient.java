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
package org.biojava.nbio.structure.align.client;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.fatcat.FatCatRigid;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.URLConnectionTools;
import org.biojava.nbio.structure.align.util.ResourceManager;
import org.biojava.nbio.structure.align.xml.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLEncoder;
import java.util.Random;
import java.util.SortedSet;
import java.util.TreeSet;

public class JFatCatClient {
	private final static Logger logger = LoggerFactory.getLogger(JFatCatClient.class);

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

	public static void main(String[] args) throws Exception {
		//System.out.println(hasPrecalculatedResult("http://source.rcsb.org/jfatcatserver/align/", "jCE Circular Permutation", "1CDG.A", "1TIM.A"));
		AtomCache cache = new AtomCache();
		String name1= "2W72.A";
		String name2= "1D2Z.D";

		Atom[] ca1 = cache.getAtoms(name1);
		Atom[] ca2 = cache.getAtoms(name2);

		int timeout = 10000;

		String testServer = "http://source.rcsb.org/jfatcatserver/align/";

		System.out.println(getAFPChainFromServer(testServer, FatCatRigid.algorithmName, name1, name2, ca1, ca2, timeout));

		PdbPairsMessage msg = getPdbPairs(testServer, 1, "test");

		System.out.println(msg);

		System.out.println(getRepresentatives(FarmJobParameters.DEFAULT_SERVER_URL, 40));
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

			InputStream stream = URLConnectionTools.getInputStream(url,timeout);

			String xml = null;

			if ( stream != null) {

				xml = convertStreamToString(stream);
				logger.info(" has PrecalcResults got XML from server: " + xml);
				HasResultXMLConverter conv = new HasResultXMLConverter();
				hasResults = conv.fromXML(xml);
			}

		} catch (IOException e){
			// log error and return false
			logger.error("error in JFatCatClient: getAFPChainFromServer",e);
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

			InputStream stream = URLConnectionTools.getInputStream(url,timeout);

			String xml = null;

			if ( stream != null) {

				xml = convertStreamToString(stream);
				//System.out.println("got XML from server: " + xml);
				PositionInQueueXMLConverter conv = new PositionInQueueXMLConverter();
				position = conv.fromXML(xml);
			}

		} catch (IOException e){
			logger.error("error in JFatCatClient: getAFPChainFromServer",e); // TODO dmyersturnbull: method should throw; we shouldn't catch here
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
			logger.info("requesting alignment from server..."  + url);
			// have a short timeout for this...
			// 5 sec
			InputStream stream = URLConnectionTools.getInputStream(url,timeout);

			String xml = null;

			if ( stream != null) {

				xml = convertStreamToString(stream);
			}
			if (xml != null) {

				return AFPChainXMLParser.fromXML (xml, name1, name2, ca1, ca2);

			} else {
				return null;
			}
			// TODO dmyersturnbull: method should throw; we shouldn't catch here
		} catch (IOException e){
			logger.error("error in JFatCatClient: getAFPChainFromServer",e);
		} catch (StructureException e) {
			logger.error("error in JFatCatClient: getAFPChainFromServer",e);
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
			logger.error("Couldn't convert stream to string", e); // TODO dmyersturnbull: method should throw; we shouldn't catch here
		} finally {
			try {
				stream.close();
			} catch (IOException e) {
				logger.warn("Can't close stream", e);
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
				InputStream response = URLConnectionTools.doPOST(url, multiXML,timeout);
				responseS = convertStreamToString(response);
				submitted = true;
				if (! responseS.contains("OK"))
					logger.error("server returned " + responseS);

				// server is busy... wait a bit and try again
				if ( responseS.startsWith(COME_BACK_LATER)){
					submitted = false;
				}

			} catch (Exception e){
				logger.error("Error in JFatCatClient: while sending results back to server",e);

				try {
					int randomSleep = getRandomSleepTime();
					logger.warn("sleeping " + (randomSleep/1000) + " sec.");
					Thread.sleep(randomSleep);
				} catch (InterruptedException ex){
					logger.warn("Interrupted while sleeping", ex);
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


	public static final void sendAFPChainToServer(String serverLocation, AFPChain afpChain,Atom[] ca1, Atom[] ca2) throws JobKillException
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

			InputStream response = URLConnectionTools.doPOST(url, xml,timeout);

			logger.debug("got response: {}", convertStreamToString(response));

			if ( xml.startsWith("KILL_JOB")){
				throw new JobKillException("Server responded with KILL message.");
			}

		} catch (IOException e){
			logger.error("error in JFatCatClient: sendAFPChainToServer",e);
		}

	}

	public static final int getTimeout(){
		String timeoutS = resourceManager.getString("connection.timeout");
		int timeout = 60000;

		try {
			timeout = Integer.parseInt(timeoutS);
		} catch (NumberFormatException ex ){
			logger.error("Bad connection.timeout parameter",ex);
		}
		return timeout;
	}


	public static final PdbPairsMessage getPdbPairs(String url,int nrPair, String username) throws IOException, JobKillException {


		String urlS= url + "getPairs?" + "nrPairs=" + nrPair + "&username=" + URLEncoder.encode(username, "UTF-8");
		int timeout = getTimeout();

		PdbPairsMessage msg = null;
		logger.info("requesting {}", urlS);
		URL serverUrl = new URL(urlS);
		// we are very tolerant with requesting a set of pairs, since we probably depend on getting it to get work started...
		// 1 min...
		InputStream stream = URLConnectionTools.getInputStream(serverUrl,timeout);
		String xml = null;

		if ( stream != null) {

			xml = convertStreamToString(stream);
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

		logger.info("Fetching representatives from "+u);
		try {
			URL url = new URL(u);

			InputStream stream = URLConnectionTools.getInputStream(url,timeout);

			String xml = null;

			if ( stream != null) {

				xml = convertStreamToString(stream);
			}
			if (xml != null) {
				representatives = RepresentativeXMLConverter.fromXML(xml);
			}
		} catch (IOException e){ // TODO dmyersturnbull: method should throw; we shouldn't catch here
			logger.error("Error fetching representatives",e);
		}
		return representatives;
	}



}
