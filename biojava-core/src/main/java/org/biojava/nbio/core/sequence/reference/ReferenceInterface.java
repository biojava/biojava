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
