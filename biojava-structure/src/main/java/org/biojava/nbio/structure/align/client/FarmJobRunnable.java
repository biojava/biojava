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
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignmentFactory;
import org.biojava.nbio.structure.align.ce.CeCPMain;
import org.biojava.nbio.structure.align.ce.CeMain;
import org.biojava.nbio.structure.align.events.AlignmentProgressListener;
import org.biojava.nbio.structure.align.fatcat.FatCatFlexible;
import org.biojava.nbio.structure.align.fatcat.FatCatRigid;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AFPChainScorer;
import org.biojava.nbio.structure.align.util.AlignmentTools;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.ResourceManager;
import org.biojava.nbio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.nbio.structure.align.xml.PdbPairsMessage;
import org.biojava.nbio.structure.domain.RemotePDPProvider;
import org.biojava.nbio.structure.io.LocalPDBDirectory.FetchBehavior;
import org.biojava.nbio.structure.scop.RemoteScopInstallation;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.biojava.nbio.core.util.FlatFileCache;
import org.biojava.nbio.core.util.PrettyXMLWriter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.util.*;




/** Contains the single thread for a job that can run multiple alignments.
 *
 * @author Andreas Prlic
 *
 */
public class FarmJobRunnable implements Runnable {

	private static final Logger logger = LoggerFactory.getLogger(FarmJobRunnable.class);


	//private static final int DEFAULT_PAIR_FETCH_DELAY   = 30000;
	//private static final String CONNECTION_PAIR_DELAY   = "connection.pair.delay";
	private static final String JFATCAT_NAME            = "jfatcat.name";
	private static final String JFATCAT_VERSION         = "jfatcat.version";

	private static ResourceManager resourceManager = ResourceManager.getResourceManager("jfatcat");


	//private static DateFormat dateFormat = new SimpleDateFormat("MMMM dd, yyyy h:mm a",Locale.US);

	FarmJobParameters params;

	String prevName1;
	Atom[] ca1 ;


	long startTime;
	long maxTime;
	int maxNrAlignments;
	int alignmentsCalculated;

	boolean waitForAlignments;

	private static final String randomUsername = getRandomUsername();

	boolean terminated ;

	List<AlignmentProgressListener> progressListeners;
	CountProgressListener counter ;

	String userName = null;
	protected AtomCache cache;

	boolean verbose = false; // TODO dmyersturnbull: we should probably remove this in favor of SLF4J
	String version = null;

	private static final String alignURL = "/align/";
	public FarmJobRunnable(FarmJobParameters params){
		terminated = false;
		this.params = params;
		verbose = false;

		// multiple farm jobs share the same SoftHashMap for caching coordinates
		cache = new AtomCache( params.getPdbFilePath(), params.getCacheFilePath());


		if ( params.getServer()!= null && (!params.getServer().equals("") ) ) {

			RemotePDPProvider pdpprovider = new RemotePDPProvider();

			String serverURL = params.getServer();
			if ( ! serverURL.endsWith("/"))
				serverURL += "/";

			if (  serverURL.endsWith(alignURL)) {
				serverURL = serverURL.substring(0,serverURL.length()-alignURL.length());
			}

			pdpprovider.setServer(serverURL+"/domains/");

			cache.setPdpprovider(pdpprovider);

			RemoteScopInstallation scop = new RemoteScopInstallation();

			scop.setServer(serverURL+"/domains/");
			ScopFactory.setScopDatabase(scop);

		}

		cache.setFetchBehavior(FetchBehavior.FETCH_FILES);

		maxNrAlignments = params.getNrAlignments();
		progressListeners = null;
		if (params.getUsername() == null) {
			userName = randomUsername;
		} else {
			userName = params.getUsername();
		}
		counter  = new CountProgressListener();
		addAlignmentProgressListener(counter);
		waitForAlignments = true;

		if ( params.isVerbose()){
			verbose = true;
		}
	}

	public void addAlignmentProgressListener(AlignmentProgressListener listener){

		if (progressListeners == null)
			progressListeners = new ArrayList<AlignmentProgressListener>();

		progressListeners.add(listener);
	}

	public void clearListeners(){
		if ( progressListeners == null)
			return;
		progressListeners.clear();
		progressListeners = null;
	}

	protected static String getRandomUsername(){
		String name = "";
		try {
			InetAddress i = InetAddress.getLocalHost();
			name += i.getHostAddress();
			name += "_";
		} catch (UnknownHostException e){
			throw new RuntimeException(e);
		}
		name +=  UUID.randomUUID();

		return name;

	}

	@Override
	public void run() {

		// Retrieve resource
		String appVersion = resourceManager.getString(JFATCAT_VERSION);
		String appName    = resourceManager.getString(JFATCAT_NAME);
		logger.info("{} version: {}", appName, appVersion);


		startTime = System.currentTimeMillis();
		// -t ime is in seconds.
		long maxSec = params.getTime();

		if ( maxSec < 5 )
			maxTime = Long.MAX_VALUE;
		else
			maxTime = startTime + params.getTime() * 1000;

		terminated = false;

		alignmentsCalculated = 0;

		maxNrAlignments = params.getNrAlignments();

		if ( maxNrAlignments < 0 ){
			maxNrAlignments = Integer.MAX_VALUE;
		}

		logger.info("running job for max {} alignments", maxNrAlignments);


		while (! terminated){

			// talk to server
			// get list of alignments to run
			// if maxNrAlignments > 100 we split up the calculations in chunks of 100.
			// otherwise we request all of them at once.
			// we request
			PdbPairsMessage msg = getAlignmentPairsFromServer();
			if ( msg == null) {
				logger.error("Got null instead of alignment pairs from server.");
				randomSleep();
				continue;
			}
			SortedSet<PdbPair> alignmentPairs = msg.getPairs();
			logger.debug("{}: Server responded with {} pairs.", userName, alignmentPairs.size());
			List<String> results = new ArrayList<String>();

			String algorithmName = msg.getMethod();
			if ( version == null) {
				setVersion(algorithmName);

			}
			for(PdbPair pair : alignmentPairs){

				if ( terminated)
					break;

				long now = System.currentTimeMillis();
				if ( now >= maxTime)  {
					terminated = true;
					break;
				}

				if ( alignmentsCalculated >= maxNrAlignments) {
					terminated = true;
					break;
				}


				String name1 = pair.getName1();
				String name2 = pair.getName2();

				if ( progressListeners != null)
					notifyStartAlignment(name1,name2);


				try {
					String resultXML = alignPair(name1, name2,algorithmName);

					if ( progressListeners != null)
						notifyEndAlignment();

					results.add(resultXML);

				} catch (Exception e){
					logger.error("Problem aligning {} with {}", name1, name2, e);

					StringWriter sw = new StringWriter();

					PrettyXMLWriter xml = new PrettyXMLWriter(new PrintWriter(sw));
					try {
						xml.openTag("AFPChain");

						xml.attribute("name1", name1);
						xml.attribute("name2", name2);
						xml.attribute("error", e.getMessage());
						xml.closeTag("AFPChain");
					} catch(IOException ex){
						logger.error("Error occured converting alignment for {} and {} to XML", name1, name2, ex);
					}

					if ( progressListeners != null)
						notifyEndAlignment();

					results.add(sw.toString());


				}

				alignmentsCalculated++;

			}

			// send results back to server
			sendResultsToServer(results);

			long end = System.currentTimeMillis();
			if ( end >= maxTime)  {
				logger.info("OK end of job: reached maxTime {}", maxTime);
				terminated = true;

			}

			if ( alignmentsCalculated >= maxNrAlignments) {
				logger.info("OK end of job: reached maxNrAlignments", maxNrAlignments);
				terminated = true;

			}

			long tdiff = (end - startTime);
			if ( tdiff != 0) {

				logger.info(userName + String.format(": job has run for :  %.2f", (tdiff) / 1000.0 / 60) + " min.");
				logger.info("{}: total nr of alignments calculated: {}", userName, alignmentsCalculated);
				if ( alignmentsCalculated > 0)
					logger.info(userName + String.format(": average time / alignment: %.2f", (tdiff / alignmentsCalculated / 1000.0)) + " sec.");
			}
		}

		logger.info(userName + ": jFATCAT job result: " + counter);

		// clean up in the end...
		clearListeners();

		cache.notifyShutdown();

	}


	private void setVersion(String algorithmName) {
		StructureAlignment algorithm;
		try {
			algorithm = StructureAlignmentFactory.getAlgorithm(algorithmName);
			version = algorithm.getVersion();
		} catch (StructureException e) {
			throw new RuntimeException("Couldn't set version for algorithm \"" + algorithmName + "\"", e);
//			version = resourceManager.getString(JFATCAT_VERSION); // dmyersturnbull: was this
		}


	}

	private void notifyStartAlignment(String name1, String name2) {
		if ( progressListeners != null){
			for (AlignmentProgressListener li : progressListeners){
				li.alignmentStarted(name1, name2);
			}
		}
	}

	private void notifyEndAlignment(){
		if ( progressListeners != null){
			for (AlignmentProgressListener li : progressListeners){
				li.alignmentEnded();

			}
		}
	}

	private void notifyRequestingAlignments(int nrAlignments){
		if ( progressListeners != null){
			for (AlignmentProgressListener li : progressListeners){
				li.requestingAlignmentsFromServer(nrAlignments);

			}
		}
	}

	private void notifySubmittingAlignments(int nrAlignments, String message){
		if ( progressListeners != null){
			for (AlignmentProgressListener li : progressListeners){
				li.sentResultsToServer(nrAlignments,message);

			}
		}
	}


	public String alignPair(String name1, String name2)
		throws StructureException, IOException {
		return alignPair(name1, name2, FatCatRigid.algorithmName);
	}

	public String alignPair(String name1, String name2, String algorithmName)
		throws StructureException, IOException {

		// 	make sure each thread has an independent instance of the algorithm object ...

		StructureAlignment algorithm = getAlgorithm(algorithmName);

		// we are running with default parameters

		if ( verbose ) {
			logger.debug("aligning {} against {}", name1, name2);
		}

		long startTime = System.currentTimeMillis();

		if ( prevName1 == null)
			initMaster(name1);

		if ( ! prevName1.equals(name1) ) {
			// we need to reload the master
			initMaster(name1);
		}

		// get a copy of the atoms, but clone them, since they will be rotated...
		Atom[] ca2 =  cache.getAtoms(name2);

		AFPChain afpChain = algorithm.align(ca1, ca2);

		afpChain.setName1(name1);
		afpChain.setName2(name2);

		try {
			// add tmScore
			double tmScore = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
			afpChain.setTMScore(tmScore);
		} catch (RuntimeException e){
			logger.error("ca1 size: {} ca2 length: {} {}  {}", ca1.length, ca2.length, afpChain.getName1(), afpChain.getName2(), e);

		}
		long endTime = System.currentTimeMillis();

		long calcTime = (endTime-startTime);
		if ( verbose ){
			boolean isCP = !AlignmentTools.isSequentialAlignment(afpChain, false);
			String msg = "finished alignment: " + name1 + " vs. " + name2 + " in " + (calcTime) / 1000.0 + " sec.";
			msg += " algo: " + algorithmName + " v:" + version + " " + afpChain;

			if ( isCP ) msg += "HAS A CIRCULAR PERMUTATION!!!";
			logger.debug(msg);
		}
		if (verbose){
			printMemory();
		}
		afpChain.setCalculationTime(calcTime);

		return AFPChainXMLConverter.toXML(afpChain, ca1, ca2);
	}




	private void printMemory() {
		int size = 1048576;
		long heapSize = Runtime.getRuntime().totalMemory() / size;

		// Get maximum size of heap in bytes. The heap cannot grow beyond this size.
		// Any attempt will result in an OutOfMemoryException.
		long heapMaxSize = Runtime.getRuntime().maxMemory() / size;

		// Get amount of free memory within the heap in bytes. This size will increase
		// after garbage collection and decrease as new objects are created.
		long heapFreeSize = Runtime.getRuntime().freeMemory() / size;
		StringBuilder msg = new StringBuilder();
		msg.append("  total: ").append(heapSize).append(" M");
		msg.append(" max: "). append(heapMaxSize).append(" M");
		msg.append(" free: ").append(heapFreeSize).append(" M");

		logger.debug(msg.toString());

	}

	private StructureAlignment getAlgorithm(String algorithmName) throws StructureException {


		StructureAlignment algorithm    = null;

		if ( algorithmName == null){

			algorithm = new FatCatRigid();

		} else if ( algorithmName.equalsIgnoreCase(FatCatRigid.algorithmName)){

				algorithm = new FatCatRigid();

		} else if ( algorithmName.equalsIgnoreCase(CeMain.algorithmName)){

			algorithm = new CeMain();

		} else if ( algorithmName.equalsIgnoreCase(CeCPMain.algorithmName)){

			algorithm = new CeCPMain();

		} else if ( algorithmName.equalsIgnoreCase(FatCatFlexible.algorithmName)){

			algorithm = new FatCatFlexible();

		} else {

			algorithm = StructureAlignmentFactory.getAlgorithm(algorithmName);

		}

		if ( algorithm == null) {

			algorithm = new FatCatRigid();

		}


		return algorithm;
	}

	private void initMaster(String name1) throws IOException, StructureException{

		ca1 = cache.getAtoms(name1);

		prevName1 = name1;

	}


	/** talk to centralized server and fetch all alignments to run.
	 *
	 * @return a list of pairs to align.
	 */
	protected PdbPairsMessage getAlignmentPairsFromServer() {


		String url = params.getServer();

		int nrPairs = params.getStepSize();

		if ( maxNrAlignments < nrPairs )
			nrPairs = maxNrAlignments;

		SortedSet<PdbPair> allPairs = new TreeSet<PdbPair>();

		PdbPairsMessage msg = null;


		try {

			if ( progressListeners != null)
				notifyRequestingAlignments(nrPairs);



			if ( ! waitForAlignments) {
				msg = JFatCatClient.getPdbPairs(url, nrPairs, userName);
				allPairs = msg.getPairs();

			} else {

				while (allPairs.isEmpty()) {
					msg = JFatCatClient.getPdbPairs(url, nrPairs, userName);
					allPairs = msg.getPairs();

					if (allPairs.isEmpty()) {
						randomSleep();
					}
				}
			}
		} catch ( JobKillException k ){

			logger.debug("Terminating job", k);
			terminate();

		} catch (Exception e) {
			logger.error("Error while requesting alignment pairs", e);
			// an error has occured sleep 30 sec.

			randomSleep();


		}

		return msg;
	}

	private void randomSleep() {
		try {

			int delay = JFatCatClient.getRandomSleepTime();
			logger.debug("sleeping {} sec.", delay/1000);
			Thread.sleep(delay);
		} catch (InterruptedException ex){
			logger.trace("InterruptedException occurred while sleeping", ex);
		}

	}

	protected void sendResultsToServer(List<String> results) {

		String serverLocation = params.getServer();

		if ( results.size() < 1)
			return;

		String fullXml = "<alignments>";

		for (String xml: results){
			fullXml +=xml;
		}
		fullXml += "</alignments>";
		String msg = "";
		try {
			msg = JFatCatClient.sendMultiAFPChainToServer(serverLocation,fullXml, userName, version );
		} catch (JobKillException e){
			logger.info("{} Got Job Kill message from server, terminating...", userName, e);
			terminate();
		}

		if ( progressListeners != null)
			notifySubmittingAlignments(results.size(), msg);
		logger.info("{}: Sent {} results to server. job status: {}", userName, results.size(), counter);
		logger.info("{}: fileCache size: {}", userName, FlatFileCache.size());
	}


	/** Send signal to terminate calculations
	 *
	 */
	public synchronized void terminate(){
		terminated = true;
	}

	public boolean isWaitForAlignments() {
		return waitForAlignments;
	}

	public void setWaitForAlignments(boolean waitForAlignments) {
		this.waitForAlignments = waitForAlignments;
	}



}
