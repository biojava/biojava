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
