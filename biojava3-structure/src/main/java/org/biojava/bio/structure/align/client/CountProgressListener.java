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
 * Created on Sep 15, 2009
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.client;

import org.biojava.bio.structure.align.events.AlignmentProgressListener;

public class CountProgressListener implements AlignmentProgressListener {

	int nrCalculated ;
	int nrSubmitted;
	
	public CountProgressListener(){
		nrCalculated = 0;
		nrSubmitted  = 0;
	}
	
	public void alignmentEnded() {
		nrCalculated++;

	}

	public void alignmentStarted(String name1, String name2) {
		// TODO Auto-generated method stub

	}

	public void downloadingStructures(String name) {
		// TODO Auto-generated method stub

	}

	public void logStatus(String message) {
		// TODO Auto-generated method stub

	}

	public void requestingAlignmentsFromServer(int nrAlignments) {
		// TODO Auto-generated method stub

	}

	public void sentResultsToServer(int nrAlignments, String serverMessage) {
		nrSubmitted+=nrAlignments;				
	}


	public String toString() {
		return "[nrCalculated=" + nrCalculated
				+ ", nrSubmitted=" + nrSubmitted + "]";
	}
	
	

}
