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
 * Created on DATE
 *
 */
package org.biojava.nbio.core.sequence;

import org.biojava.nbio.core.util.Equals;
import org.biojava.nbio.core.util.Hashcoder;

/**
 * Used in Sequences as the unique indentifier. If possible, set the {@link DataSource} to know the
 * source of the id. This allows a SequenceProxy to gather features or related sequences
 * Protein->Gene as an example. When parsing a Blast file it is also possible
 * to identify the type of ID
 *
 * @author Scooter Willis
 * @author Jacek Grzebyta
 */
public class AccessionID {

	private String id = null;
	private DataSource source = DataSource.LOCAL;
	private Integer version;
	private String identifier = null;

	/**
	 *
	 */

	public AccessionID(){
		id = "";

	}

	/**
	 *
	 * @param id
	 */
	public AccessionID(String id) {
		this.id = id.trim();
		this.source = DataSource.LOCAL;
	}


	/**
	 *
	 * @param id
	 * @param source
	 */
	public AccessionID(String id, DataSource source) {
		this.id = id.trim();
		this.source = source;
	}

	public AccessionID(String id, DataSource source, Integer version, String identifier) {
		this.id = id;
		this.source = source;
		this.version = version;
		this.identifier = identifier;
	}

	/**
	 * @return the id
	 */
	public String getID() {
		return id;
	}

	/**
	 * @return the source
	 */
	public DataSource getDataSource() {
		return source;
	}

	@Override
	public boolean equals(Object o) {
		boolean equals = false;
		if (Equals.classEqual(this, o)) {
			AccessionID l = (AccessionID) o;
			equals = (Equals.equal(getID(), l.getID())
					&& Equals.equal(getDataSource(), l.getDataSource())
					&& Equals.equal(getIdentifier(), l.getIdentifier())
					&& Equals.equal(getVersion(), l.getVersion()));
	}
		return equals;
	}

	@Override
	public int hashCode() {
		int r = Hashcoder.SEED;
		r = Hashcoder.hash(r, getID());
		r = Hashcoder.hash(r, getDataSource());
		r = Hashcoder.hash(r, getIdentifier());
		r = Hashcoder.hash(r, getVersion());
		return r;
	}

//   public void setDataSource(DataSource dataSource){
//       source = dataSource;
//   }


	/**
	 * In case if the {@link #getID() } is not unique keeps the id version.
	 * @return the version
	 */
	public Integer getVersion() {
		return version;
	}

	public void setVersion(Integer version) {
		this.version = version;
	}

	/**
	 * In case if {@link #getID() } in not unique keeps the alternative id, eg. NCBI GI number.
	 * 
	 * This may null.
	 *
	 * @return
	 */
	public String getIdentifier() {
		return identifier;
	}

	public void setIdentifier(String identifier) {
		this.identifier = identifier;
	}


	@Override
	public String toString() {
		return id;
	}
}
