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
public abstract class AbstractReference implements ReferenceInterface {

    /**
     * The title that retrieved from the Reference section.
     */
    private String title;

    /**
     * The authors are a list of Inventors that retrieved from the Reference section.
     */
    private String authors;

    /**
     * The journal usually contains the Publication Number, Publication Date and Assignee
     */
    private String journal;

    /**
     * The title that retrieved from the Reference section.
     *
     * @return
     */
    @Override
    public String getTitle() {
        return title;
    }

    /**
     * Set The title that retrieved from the Reference section.
     *
     * @param title
     */
    @Override
    public void setTitle(String title) {
        this.title = title;
    }

    /**
     * The authors are a list of Inventors that retrieved from the Reference section.
     *
     * @return
     */
    @Override
    public String getAuthors() {
        return authors;
    }

    /**
     * Set The authors are a list of Inventors that retrieved from the Reference section.
     *
     * @param authors
     */
    @Override
    public void setAuthors(String authors) {
        this.authors = authors;
    }

    /**
     * The journal usually contains the Publication Number, Publication Date and Assignee
     *
     * @return
     */
    @Override
    public String getJournal() {
        return journal;
    }

    /**
     * Set The journal usually contains the Publication Number, Publication Date and Assignee
     *
     * @param journal
     */
    @Override
    public void setJournal(String journal) {
        this.journal = journal;
    }
}
