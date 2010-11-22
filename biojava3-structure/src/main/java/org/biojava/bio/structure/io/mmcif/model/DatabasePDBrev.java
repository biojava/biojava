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
 * created at Apr 27, 2008
 */
package org.biojava.bio.structure.io.mmcif.model;

public class DatabasePDBrev {
	String date;
	String date_original;
	String status;
	String replaces;
	String mod_type;
	String num;

	public String toString(){
		StringBuffer buf = new StringBuffer();
		buf.append("DatabasePDBrev ");
		buf.append("mod_type :");
		buf.append(mod_type);
		buf.append(" ");
		buf.append(this.getDate());
		buf.append( " ");
		buf.append( this.getDate_original());

		return buf.toString();
	}
	public String getNum() {
		return num;
	}
	public void setNum(String num) {
		this.num = num;
	}
	public String getDate() {
		return date;
	}
	public void setDate(String date) {
		this.date = date;
	}
	public String getDate_original() {
		return date_original;
	}
	public void setDate_original(String date_original) {
		this.date_original = date_original;
	}
	public String getStatus() {
		return status;
	}
	public void setStatus(String status) {
		this.status = status;
	}
	public String getReplaces() {
		return replaces;
	}
	public void setReplaces(String replaces) {
		this.replaces = replaces;
	}
	public String getMod_type() {
		return mod_type;
	}
	public void setMod_type(String mod_type) {
		this.mod_type = mod_type;
	}


}
