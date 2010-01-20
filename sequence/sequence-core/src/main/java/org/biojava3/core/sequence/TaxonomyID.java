/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence;

/**
 *
 * @author Scooter
 */
public class TaxonomyID {

    public enum Source {

        NCBI, LOCAL
    }
    private String id = null;
    private Source source = Source.LOCAL;

    public TaxonomyID(String id, Source source) {
        this.id = id;
        this.source = source;
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
    public Source getSource() {
        return source;
    }
}
