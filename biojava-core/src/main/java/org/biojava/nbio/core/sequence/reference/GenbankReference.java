package org.biojava.nbio.core.sequence.reference;

/**
 * For Genbank format file only.
 *
 * @Author Jim Tang
 */
public class GenbankReference extends AbstractReference {

    /**
     * The authors are a list of Inventors that retrieved from the Reference section.
     */
    private String authors;

    /**
     * The title that retrieved from the Reference section.
     */
    private String title;

    /**
     * The journal usually contains the Publication Number, Publication Date and Assignee
     */
    private String journal;

    /**
     * @return
     */
    public String getAuthors() {
        return authors;
    }

    /**
     * @param authors
     */
    public void setAuthors(String authors) {
        this.authors = authors;
    }

    /**
     * @return
     */
    public String getTitle() {
        return title;
    }

    /**
     * @param title
     */
    public void setTitle(String title) {
        this.title = title;
    }

    /**
     * @return
     */
    public String getJournal() {
        return journal;
    }

    /**
     * @param journal
     */
    public void setJournal(String journal) {
        this.journal = journal;
    }
}
