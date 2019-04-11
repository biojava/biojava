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
package org.biojava.nbio.core.sequence.reference;

/**
 * @since 5.0.0
 * @Author Jim Tang
 */
public interface ReferenceInterface {

	/**
	 * Set the title that retrieved from Reference section.
	 *
	 * @param title
	 */
	void setTitle(String title);

	/**
	 * Get the title that retrieved from Reference section.
	 *
	 * @return
	 */
	String getTitle();

	/**
	 * Set the authors that retrieved from Reference section.
	 *
	 * @param authors
	 */
	void setAuthors(String authors);

	/**
	 * Get the authors that retrieved from Reference section.
	 *
	 * @return
	 */
	String getAuthors();

	/**
	 * Set the journal that retrieved from Reference section.
	 *
	 * @param journal
	 */
	void setJournal(String journal);

	/**
	 * Get the journal that retrieved from Reference section.
	 *
	 * @return
	 */
	String getJournal();

}
