package org.biojava.nbio.core.sequence.reference;

/**
 * @Author Jim Tang
 */
public abstract class AbstractReference implements ReferenceInterface {

    private String title;

    private String authors;

    private String journal;

    @Override
    public String getTitle() {
        return title;
    }

    @Override
    public void setTitle(String title) {
        this.title = title;
    }

    @Override
    public String getAuthors() {
        return authors;
    }

    @Override
    public void setAuthors(String authors) {
        this.authors = authors;
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
