
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
 * Created by andreas on 10/12/15.
 */

package org.biojava.nbio.structure.io.mmcif.model;

import java.io.Serializable;

public class DatabasePdbrevRecord implements Serializable {


	private static final long serialVersionUID = 1L;

	String rev_num;
	String type;
	String details;

	public String getRev_num() {
		return rev_num;
	}

	public void setRev_num(String rev_num) {
		this.rev_num = rev_num;
	}

	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}

	public String getDetails() {
		return details;
	}

	public void setDetails(String details) {
		this.details = details;
	}

	@Override
	public String toString() {
		return "DatabasePdbrevRecord{" +
				"rev_num='" + rev_num + '\'' +
				", type='" + type + '\'' +
				", details='" + details + '\'' +
				'}';
	}
}
