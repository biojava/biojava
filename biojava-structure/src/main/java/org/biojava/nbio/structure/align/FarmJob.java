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
package org.biojava.nbio.structure.align;

import org.biojava.nbio.structure.align.client.FarmJobParameters;
import org.biojava.nbio.structure.align.client.FarmJobRunnable;
import org.biojava.nbio.structure.align.events.AlignmentProgressListener;
import org.biojava.nbio.structure.align.util.CliTools;
import org.biojava.nbio.structure.align.util.ConfigurationException;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.scop.CachedRemoteScopInstallation;
import org.biojava.nbio.structure.scop.ScopDatabase;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.biojava.nbio.core.util.InputStreamProvider;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


/** A job as it can be run on the farm.
 *
 * @author Andreas Prlic
 *
 * for arguments see the printHelp() method.
 *
 *
 *
 */
public class FarmJob implements Runnable {

	private final static Logger logger = LoggerFactory.getLogger(FarmJob.class);

	private static final String[] mandParams = new String[] {"pdbFilePath"};

	private static final List<String> mandatoryArgs= Arrays.asList(mandParams);

	List<AlignmentProgressListener> progressListeners;
	List<FarmJobRunnable> jobs;

	FarmJobParameters params ;

	public FarmJob(){
		progressListeners = null;

		// send a flag to the PDb file loader to cache the gzip compressed files.
		System.setProperty(InputStreamProvider.CACHE_PROPERTY, "true");


	}

	public FarmJobParameters getParams() {
		return params;
	}

	public void setParams(FarmJobParameters params) {
		this.params = params;
	}

	public void addAlignmentProgressListener(AlignmentProgressListener listener){
		if (progressListeners == null)
			progressListeners = new ArrayList<AlignmentProgressListener>();

		progressListeners.add(listener);
	}

	public void clearListeners(){
		progressListeners.clear();
		progressListeners = null;
	}

	public static void main(String[] argv){

		FarmJob job = new FarmJob();

		if (argv.length  == 0 ) {
			job.printHelp();
			return;
		}

		if ( argv.length == 1){
			if (argv[0].equalsIgnoreCase("-h") || argv[0].equalsIgnoreCase("-help")|| argv[0].equalsIgnoreCase("--help")){
				job.printHelp();
				return;
			}
		}

		FarmJobParameters params = new FarmJobParameters();

		for (int i = 0 ; i < argv.length; i++){
			String arg   = argv[i];

			String value = null;
			if ( i < argv.length -1)
				value = argv[i+1];

			// if value starts with - then the arg does not have a value.
			if (value != null && value.startsWith("-"))
				value = null;
			else
				i++;


			String[] tmp = {arg,value};

			try {

				CliTools.configureBean(params, tmp);

			} catch (ConfigurationException e){

				logger.error("Exception", e);

				if ( mandatoryArgs.contains(arg) ) {
					// there must not be a ConfigurationException with mandatory arguments.
					return;
				} else {
					// but there can be with optional ...
				}
			}
		}


		if (( params.getNrAlignments() == -1) && (params.getTime() == -1)){
			logger.error("Please provide either the -time or the -nrAlignments argument!");
			return;
		}


		logger.info("Using parameters: {}", params);

		job.setParams(params);
		job.run();

	}

	@Override
	public void run(){


		// set the system wide PDB path

		String path = params.getPdbFilePath();
		System.setProperty(UserConfiguration.PDB_DIR,path);

		String cachePath = params.getCacheFilePath();
		if ( cachePath != null && ! cachePath.equals(""))
			System.setProperty(UserConfiguration.PDB_CACHE_DIR,cachePath);
		else {
			// if not provided, we use pdbFilePath as the default CACHE path
			System.setProperty(UserConfiguration.PDB_CACHE_DIR,path);
		}
		// declare SCOP to be locally cached, but fetching new stuff from remote
		ScopDatabase scop = null;
		try {
			scop = new CachedRemoteScopInstallation(true);
		} catch (IOException e) {
			throw new RuntimeException("Could not load " + CachedRemoteScopInstallation.class.getName(), e);
		}
		ScopFactory.setScopDatabase(scop);

		String username = params.getUsername();
		jobs = new ArrayList<FarmJobRunnable>();
		for ( int i = 0 ; i < params.getThreads();i++){
			logger.info("starting thread #{}", (i+1));
			FarmJobRunnable runner = new FarmJobRunnable(params);
			params.setUsername(username+"_thread_" + (i+1));
			jobs.add(runner);

			if ( progressListeners != null) {
				for (AlignmentProgressListener li : progressListeners){
					runner.addAlignmentProgressListener(li);
				}
			}


			Thread t = new Thread(runner);
			if ( ( (params.getThreads() > 1 ) && ( i < params.getThreads() - 1) )|| ( params.isRunBackground())) {
				logger.info("starting thread #{} in background.", (i + 1));
				t.start();
			} else {
				// single CPU jobs are run in the main thread and the last job is also run in the main thread
				logger.info("starting thread #{} in main thread.", (i + 1));
				t.run();
			}
		}
	}

	public void terminate(){

		logger.info("terminating jobs");

		if ( jobs == null)
			return;

		int js = jobs.size();
		logger.info("number of jobs: {}", js);


		for (FarmJobRunnable runner : jobs){
			// runner.terminate() is already synchronized
			runner.terminate();
		}

		clearListeners();
	}

	public void printHelp(){
		System.out.println("-------------------");
		System.out.println("FarmJob help:");
		System.out.println("-------------------");

		System.out.println("FarmJob accepts the following parameters:");
		System.out.println("");
		System.out.println(" Mandatory:");
		System.out.println("   -pdbFilePath (mandatory) Path to the directory in your file system that contains the PDB files.");

		System.out.println("   provide either -time or -nrAlignments. If both are provided the job stops as soon as any of the criteria has been reached.");
		System.out.println("   -time maximum number of time to run (in seconds). -1 means no time limit, but run -nrAlignment arguments. Default: " + FarmJobParameters.DEFAULT_JOB_TIME );
		System.out.println("   -nrAlignments number of alignments to calculate. Default: " + FarmJobParameters.DEFAULT_NR_ALIGNMENTS) ;
		System.out.println("");
		System.out.println(" Optional: ");
		System.out.println("   -threads number of parallel threads to calculate alignments. Should be nr. of available CPUs. Default: " + FarmJobParameters.DEFAULT_NR_THREADS);
		System.out.println("   -server the location of the server URL to talk to. Default : " + FarmJobParameters.DEFAULT_SERVER_URL);
		System.out.println("   -username a unique name that can be given to this client. Can be used to give credit for who is doing the calculations. Default: IP and a random id");
		System.out.println("   -stepSize the number of pairs to be requsted from server. Default: " + FarmJobParameters.DEFAULT_BATCH_SIZE);
	}
}
