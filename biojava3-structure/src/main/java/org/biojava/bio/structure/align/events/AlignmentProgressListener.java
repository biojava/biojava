package org.biojava.bio.structure.align.events;

public interface AlignmentProgressListener {

	public void alignmentStarted(String name1, String name2);
	
	public void alignmentEnded();
	
	public void logStatus(String message);
	
	public void downloadingStructures(String name);
	
	public void requestingAlignmentsFromServer(int nrAlignments);
	
	public void sentResultsToServer(int nrAlignments,String serverMessage);
	
	
}
