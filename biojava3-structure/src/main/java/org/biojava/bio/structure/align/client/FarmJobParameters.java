package org.biojava.bio.structure.align.client;

import org.biojava.bio.structure.align.util.ResourceManager;

public class FarmJobParameters {

	
	public static final int DEFAULT_JOB_TIME = -1;
	public static final int DEFAULT_NR_ALIGNMENTS = -1;
	public static final int DEFAULT_NR_THREADS = 1;
	public static final String DEFAULT_SERVER_URL;
	private static ResourceManager resourceManager;
	static {
		resourceManager = ResourceManager.getResourceManager("jfatcat");
		String server = resourceManager.getString("server.url");
		DEFAULT_SERVER_URL = server;
	}
	public static final String DEFAULT_PDB_PATH = "/tmp/";
	public static final boolean DEFAULT_DIR_SPLIT = true;
	public static final int DEFAULT_BATCH_SIZE         = 100;
	
	private static final String DEFAULT_BATCH_SIZE_PROP = "request.pair.size";
	
	int nrAlignments;	
	int time;
	int threads;
	String server;
	String pdbFilePath;
	boolean pdbDirSplit;
	String username;
	boolean runBackground;
	boolean verbose;
	boolean updateRemediatedFiles;
	int stepSize;
	String cacheFilePath;

		
	public FarmJobParameters(){
		nrAlignments = DEFAULT_NR_ALIGNMENTS;
		time = DEFAULT_JOB_TIME;
		threads = DEFAULT_NR_THREADS;
		server = DEFAULT_SERVER_URL;
		pdbFilePath = DEFAULT_PDB_PATH;
		pdbDirSplit = DEFAULT_DIR_SPLIT;
		runBackground = false;
		cacheFilePath = DEFAULT_PDB_PATH;
		updateRemediatedFiles = false;
		String nrPairsProp = resourceManager.getString(DEFAULT_BATCH_SIZE_PROP);
		
		stepSize = 	DEFAULT_BATCH_SIZE;
		
		username = FarmJobRunnable.getRandomUsername();
		if ( nrPairsProp != null){
			try {
				stepSize = Integer.parseInt(nrPairsProp);
			} catch (NumberFormatException ex){
				ex.printStackTrace();
			}
		}
		
	}
	
	public String getPdbFilePath() {
		return pdbFilePath;
	}

	public void setPdbFilePath(String pdbFilePath) {
		this.pdbFilePath = pdbFilePath;
	}
	public String getCacheFilePath() {
		return cacheFilePath;
	}

	public void setCacheFilePath(String cacheFilePath) {
		this.cacheFilePath = cacheFilePath;
	}

	public int getNrAlignments() {
		return nrAlignments;
	}


	public void setNrAlignments(int nrAlignments) {
		this.nrAlignments = nrAlignments;
	}


	public int getTime() {
		return time;
	}

	public void setTime(int time) {
		this.time = time;
	}

	public int getThreads() {
		return threads;
	}

	public void setThreads(int threads) {
		this.threads = threads;
	}

	public String getServer() {
		return server;
	}

	public void setServer(String server) {
		this.server = server;
	}
	
	public boolean isPdbDirSplit() {
		return pdbDirSplit;
	}

	public void setPdbDirSplit(boolean pdbDirSplit) {
		this.pdbDirSplit = pdbDirSplit;
	}
	
	public String getUsername() {
		return username;
	}
	public void setUsername(String username) {
		this.username = username;
	}
	
	/** Flag if a job that only runs one parallell job should be run in its own thread or in the main thread.
	 * For User interface related apps should be set to true. Default: false;
	 * @return flag
	 */
	public boolean isRunBackground() {
		return runBackground;
	}
	public void setRunBackground(boolean runBackground) {
		this.runBackground = runBackground;
	}


	/** how many pairs should be requested for alignment from server?
	 * 
	 * @return stepsize
	 */
	public int getStepSize() {
		return stepSize;
	}

	public void setStepSize(int stepSize) {
		this.stepSize = stepSize;
	}

	
	/** Flag if the job should be run in verbose mode. Default: false
	 * 
	 * @return flag if the job should be run in verbose mode
	 */
	public boolean isVerbose() {
		return verbose;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	public boolean isUpdateRemediatedFiles() {
		return updateRemediatedFiles;
	}

	public void setUpdateRemediatedFiles(boolean updateRemediatedFiles) {
		this.updateRemediatedFiles = updateRemediatedFiles;
	}

	@Override
	public String toString() {
		return "FarmJobParameters [nrAlignments=" + nrAlignments + ", time="
				+ time + ", threads=" + threads + ", server=" + server
				+ ", pdbFilePath=" + pdbFilePath + ", pdbDirSplit="
				+ pdbDirSplit + ", username=" + username + ", runBackground="
				+ runBackground + ", verbose=" + verbose
				+ ", updateRemediatedFiles=" + updateRemediatedFiles
				+ ", stepSize=" + stepSize + ", cacheFilePath=" + cacheFilePath
				+ "]";
	}
	
	
	
}
