package org.biojava.bio.structure.align.client;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.net.InetAddress;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Locale;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.UUID;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.events.AlignmentProgressListener;
import org.biojava.bio.structure.align.fatcat.FatCatFlexible;
import org.biojava.bio.structure.align.fatcat.FatCatRigid;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AFPChainScorer;
import org.biojava.bio.structure.align.util.AlignmentTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.ResourceManager;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.bio.structure.align.xml.PdbPairsMessage;
import org.biojava.bio.structure.domain.RemotePDPProvider;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.scop.RemoteScopInstallation;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.core.util.FlatFileCache;
import org.biojava3.core.util.PrettyXMLWriter;




/** Contains the single thread for a job that can run multiple alignments.
 * 
 * @author Andreas Prlic
 *
 */
public class FarmJobRunnable implements Runnable {


	//private static final int DEFAULT_PAIR_FETCH_DELAY   = 30000;
	//private static final String CONNECTION_PAIR_DELAY   = "connection.pair.delay";
	private static final String JFATCAT_NAME            = "jfatcat.name";
	private static final String JFATCAT_VERSION         = "jfatcat.version";

	private static ResourceManager resourceManager = ResourceManager.getResourceManager("jfatcat");


	private static DateFormat dateFormat = new SimpleDateFormat("MMMM dd, yyyy h:mm a",Locale.US);

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

	boolean verbose = false;
	String version = null;
	
	private static final String alignURL = "/align/";
	public FarmJobRunnable(FarmJobParameters params){
		terminated = false;
		this.params = params;
		verbose = false;

		// multiple farm jobs share the same SoftHashMap for caching coordinates
		cache = new AtomCache( params.getPdbFilePath(), params.isPdbDirSplit());
		
			
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
		
		
		// enforce to replace remediated files with new versions...
		FileParsingParameters fparams = cache.getFileParsingParams();
		fparams.setUpdateRemediatedFiles(params.isUpdateRemediatedFiles());
		
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
		} catch (Exception e){
			//e.printStackTrace();
		}
		name +=  UUID.randomUUID();

		return name;

	}

	public static void log(String message){
		StringBuffer buf = new StringBuffer();

		buf.append("[");
		Date date = new Date();
		buf.append(dateFormat.format(date));
		buf.append("] ");
		buf.append(message);
		System.out.println(buf.toString());
	}

	public void run() {

		// Retrieve resource
		String appVersion = resourceManager.getString(JFATCAT_VERSION);
		String appName    = resourceManager.getString(JFATCAT_NAME);
		log(appName + " version:" + appVersion);


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

		log("running job for max: " + maxNrAlignments + " alignments");


		while (! terminated){

			// talk to server
			// get list of alignments to run
			// if maxNrAlignments > 100 we split up the calculations in junks of 100.
			// otherwise we request all of them at once.
			// we request
			PdbPairsMessage msg = getAlignmentPairsFromServer();
			if ( msg == null) {
				System.err.println("Got null instead of alignment pairs from server.");
				randomSleep();
				continue;
			}
			SortedSet<PdbPair> alignmentPairs = msg.getPairs(); 
			log(userName+": Server responded with " + alignmentPairs.size() + " pairs.");
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
					//System.out.println("calculating alignent: " + name1 + "  " + name2);
					String resultXML = alignPair(name1, name2,algorithmName);

					if ( progressListeners != null)
						notifyEndAlignment();

					//System.out.println("got XML: " + resultXML);
					results.add(resultXML);

				} catch (Exception e){
					//if (e.getMessage() == null)
					System.err.println("Problem aligning " + name1 + " " + name2);
					e.printStackTrace();
					// log that an exception has occurred and send it back to server!1
					log("Error: " + e.getMessage() + " while aligning " + name1 + " vs. " + name2);
					System.err.println(e.getMessage());
					//e.printStackTrace();

					StringWriter sw = new StringWriter();
					PrintWriter writer = new PrintWriter(sw);

					PrettyXMLWriter xml = new PrettyXMLWriter(writer);
					try {
						xml.openTag("AFPChain");

						xml.attribute("name1", name1);
						xml.attribute("name2", name2);
						xml.attribute("error", e.getMessage());
						xml.closeTag("AFPChain");
					} catch(IOException ex){
						ex.printStackTrace();
					}

					if ( progressListeners != null)
						notifyEndAlignment();

					results.add(sw.toString());


				}

				alignmentsCalculated++;

			}

			// send results back to server
			sendResultsToServer(results);
			//log("sent results to server: " + counter.toString());

			long end = System.currentTimeMillis();
			if ( end >= maxTime)  {
				System.out.println("OK end of job: reached maxTime.");
				terminated = true;

			}

			if ( alignmentsCalculated >= maxNrAlignments) {
				System.out.println("OK end of job: reached maxNrAlignments");
				terminated = true;

			}		

			long tdiff = (end - startTime);
			if ( tdiff != 0) {

				log(userName + String.format(": job has run for :  %.2f", ( tdiff)/1000.0/60) + " min.");
				log(userName + ": total nr of alignments calculated: " +alignmentsCalculated );
				if ( alignmentsCalculated > 0)
					log(userName + String.format(": average time / alignment: %.2f", ( tdiff / alignmentsCalculated / 1000.0 )) + " sec.");
			}
		}	

		log(userName+": jFATCAT job result: " + counter.toString());

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
			// TODO Auto-generated catch block
			e.printStackTrace();
			version = resourceManager.getString(JFATCAT_VERSION);
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
		//StructureAlignment fatCatRigid    = new FatCatRigid();

		// we are running with default parameters

		if ( verbose ) {
			log("aligning " + name1 + " vs. " + name2);
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

		//		if ( FatCatAligner.printTimeStamps){
		


		//			System.out.println("time to get " + name1 + " " + name2 + " : " + (endTime-startTime) / 1000.0 + " sec.");
		//		}
		AFPChain afpChain = algorithm.align(ca1, ca2);

		afpChain.setName1(name1);
		afpChain.setName2(name2);

		try {
			// add tmScore
			double tmScore = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
			afpChain.setTMScore(tmScore);
		} catch (Exception e){
			e.printStackTrace();
			System.out.println("ca1 size:" + ca1.length + " ca2 length: " + ca2.length + " " + afpChain.getName1() + " " + afpChain.getName2());
			
		}
		long endTime = System.currentTimeMillis();
		
		long calcTime = (endTime-startTime);
		if ( verbose ){
			boolean isCP = !AlignmentTools.isSequentialAlignment(afpChain, false);
			String msg = "finished alignment: " + name1 + " vs. " + name2 + " in " + (calcTime) / 1000.0 + " sec.";
			msg += " algo: " + algorithmName + " v:" + version + " " + afpChain;
			
			if ( isCP ) msg += "HAS A CIRCULAR PERMUTATION!!!";
			log(msg);
		}
		if (verbose){
			printMemory();
		}
		afpChain.setCalculationTime(calcTime);

		String xml = AFPChainXMLConverter.toXML(afpChain, ca1, ca2);
		return xml;
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
		StringBuffer msg = new StringBuffer();
		msg.append("  total: ").append(heapSize).append(" M");
		msg.append(" max: "). append(heapMaxSize).append(" M");
		msg.append(" free: ").append(heapFreeSize).append(" M");
		
		System.out.println(msg.toString());
		
	}

	private StructureAlignment getAlgorithm(String algorithmName) {
		
	
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
			
			try {
				
				algorithm = StructureAlignmentFactory.getAlgorithm(algorithmName);
		
			} catch (StructureException ex){
				ex.printStackTrace();
			}
		}
		
		if ( algorithm == null) {
			
			algorithm = new FatCatRigid();
			
		}
		
		
		return algorithm;
	}

	private void initMaster(String name1) throws IOException, StructureException{
		//AtomCache cache = AtomCache.getInstance();

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

				while (allPairs.size() == 0) {
					msg = JFatCatClient.getPdbPairs(url, nrPairs, userName);
					allPairs = msg.getPairs();

					if ( allPairs.size() == 0 ) {
						randomSleep();
					}
				}
			}
		} catch ( JobKillException k ){

			terminate();

		} catch (Exception e) {
			if ( e.getMessage() == null)
				e.printStackTrace();
			System.err.println("Error while requesting alignment pairs: " + e.getMessage());
			// an error has occured sleep 30 sec.

			randomSleep();


		}

		return msg;
	}

	private void randomSleep() {
		try {

			int delay = JFatCatClient.getRandomSleepTime();
			System.err.println("sleeping "+ delay/1000 + " sec.");
			Thread.sleep(delay);
		} catch (InterruptedException ex){
			ex.printStackTrace();
		}
		
	}

	protected void sendResultsToServer(List<String> results) {

		String serverLocation = params.getServer();

		if ( results.size() < 1)
			return;

		//System.out.println("sending " + results.size() + " results back to server");

		String fullXml = "<alignments>";

		for (String xml: results){
			fullXml +=xml;
		}
		fullXml += "</alignments>";
		String msg = "";
		try {
			msg = JFatCatClient.sendMultiAFPChainToServer(serverLocation,fullXml, userName, version );
		} catch (JobKillException e){
			log(userName+ " Got Job Kill message from server, terminating...");
			e.printStackTrace();
			terminate();
		}

		if ( progressListeners != null)
			notifySubmittingAlignments(results.size(), msg);
		log (userName + ": Sent " + results.size() +" results to server. job status:" + counter.toString());
		log (userName + ": fileCache size:" + FlatFileCache.getInstance().size());
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
