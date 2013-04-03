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
 * Created on 2011-11-20
 *
 */

package org.biojava3.ws.alignment.qblast;

/**
 * Information about QBlast search job
 * 
 * @author Gediminas Rimsa
 */
public class BlastJob {
	private String id;
	private long startTimestamp;
	private long expectedExecutionTime;

	/**
	 * Request id (RID) as received from QBlast server
	 * @return
	 */
	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public long getStartTimestamp() {
		return startTimestamp;
	}

	public void setStartTimestamp(long startTimestamp) {
		this.startTimestamp = startTimestamp;
	}

	public long getExpectedExecutionTime() {
		return expectedExecutionTime;
	}

	public void setExpectedExecutionTime(long expectedExecutionTime) {
		this.expectedExecutionTime = expectedExecutionTime;
	}

	@Override
	public String toString() {
		StringBuilder builder = new StringBuilder();
		builder.append("BlastJob [id=");
		builder.append(id);
		builder.append(", startTimestamp=");
		builder.append(startTimestamp);
		builder.append(", expectedExecutionTime=");
		builder.append(expectedExecutionTime);
		builder.append("]");
		return builder.toString();
	}

}
