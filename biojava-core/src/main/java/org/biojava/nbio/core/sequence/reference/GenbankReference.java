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
 * For Genbank format file only.
 *
 * @since 5.0.0
 * @Author Jim Tang
 */
public class GenbankReference extends AbstractReference {

	private String authors;

	private String title;

	private String journal;

	@Override
	public String getAuthors() {
		return authors;
	}

	@Override
	public void setAuthors(String authors) {
		this.authors = authors;
	}

	@Override
	public String getTitle() {
		return title;
	}

	@Override
	public void setTitle(String title) {
		this.title = title;
	}

	@Override
	public String getJournal() {
		return journal;
	}

	@Override
	public void setJournal(String journal) {
		this.journal = journal;
	}
}
