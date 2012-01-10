/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava.bio.structure;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * PDB-specific
 * @author Jules Jacobsen
 * @since 1.7
 */
public class JournalArticle implements Serializable{

    /**
    *
    */
   private static final long serialVersionUID = 5062668226159515468L;
   private List<Author> authorList = new ArrayList<Author>();
    private List<Author> editorList = new ArrayList<Author>();
    private String title = "";
    private String ref = "";
    private String journalName = "";
    private String volume;
    private String startPage;
    private int publicationDate;
    private String publisher = "";
    private String refn = "";
    private String pmid = "";
    private String doi = "";
    private boolean published = false;

    /**
     * Get the list of Authors of the JournalArticle
     *
     * @return the List of authors
     */
    public List<Author> getAuthorList() {
        return authorList;
    }

    public void setAuthorList(List<Author> authorList) {
        this.authorList = authorList;
    }

    /**
     * Get the list of editors of the JournalArticle
     *
     * @return the List of editors
     */
    public List<Author> getEditorList() {
        return editorList;
    }

    public void setEditorList(List<Author> editorList) {
        this.editorList = editorList;
    }

    /**
     * Get the value of DOI field.
     * For example: 10.1016/S0969-2126(02)00720-7
     *
     * @return the value of doi
     */
    public String getDoi() {
        return doi;
    }

    /**
     * Set the value of doi
     *
     * @param doi new value of doi
     */
    public void setDoi(String doi) {
        this.doi = doi;
    }

    /**
     * Sets the publication state of a JournalArticle - TO BE PUBLISHED == false
     * @param state
     */
    public void setIsPublished(Boolean state) {
        this.published = state;
    }
    /**
     * Get the value of PMID field.
     * For example: 12005435
     *
     * @return the value of pmid
     */
    public String getPmid() {
        return pmid;
    }

    /**
     * Set the value of pmid
     *
     * @param pmid new value of pmid
     */
    public void setPmid(String pmid) {
        this.pmid = pmid;
    }

    /**
     * Get the value of REF field.
     * For example: TO BE PUBLISHED
     *
     * @return the value of ref
     */
    public String getRef() {
        return ref;
    }

    /**
     * Set the value of the ref.
     *
     * @param ref new value of ref
     */
    public void setRef(String ref) {
        this.ref = ref;
    }

    /**
     * Get the value of REFN field.
     * For example: ISSN 0969-2126
     *
     * @return the value of ref
     */
    public String getRefn() {
        return refn;
    }

    /**
     * Set the value of the refn
     *
     * @param refn new value of refn
     */
    public void setRefn(String refn) {
        this.refn = refn;
    }

    /**
     * Get the value of title
     *
     * @return the value of title
     */
    public String getTitle() {
        return title;
    }

    /**
     * Set the value of title
     *
     * @param title new value of title
     */
    public void setTitle(String title) {
        this.title = title;
    }

    public String getJournalName() {
        return journalName;
    }

    public void setJournalName(String journalName) {
        this.journalName = journalName;
    }

    public int getPublicationDate() {
        return publicationDate;
    }

    public void setPublicationDate(int publicationDate) {
        this.publicationDate = publicationDate;
    }

    public boolean isPublished() {
        return published;
    }

    public void setPublished(boolean published) {
        this.published = published;
    }

    public String getPublisher() {
        return publisher;
    }

    public void setPublisher(String publisher) {
        this.publisher = publisher;
    }

    public String getStartPage() {
        return startPage;
    }

    public void setStartPage(String startPage) {
        this.startPage = startPage;
    }

    public String getVolume() {
        return volume;
    }

    public void setVolume(String volume) {
        this.volume = volume;
    }

    @Override
    public String toString() {
//        JRNL        AUTH   M.HAMMEL,G.SFYROERA,D.RICKLIN,P.MAGOTTI,
//        JRNL        AUTH 2 J.D.LAMBRIS,B.V.GEISBRECHT
//        JRNL        TITL   A STRUCTURAL BASIS FOR COMPLEMENT INHIBITION BY
//        JRNL        TITL 2 STAPHYLOCOCCUS AUREUS.
//        JRNL        REF    NAT.IMMUNOL.                  V.   8   430 2007
//        JRNL        REFN                   ISSN 1529-2908
//        JRNL        PMID   17351618
//        JRNL        DOI    10.1038/NI1450
        String eol = System.getProperty("line.separator");

        StringBuffer jrnlString = new StringBuffer();

        StringBuffer authString = new StringBuffer("JRNL        AUTH   ");
        StringBuffer titlString = new StringBuffer("JRNL        TITL   ");
        StringBuffer editString = new StringBuffer("JRNL        EDIT   ");
        StringBuffer refString = new StringBuffer("JRNL        REF    ");
        StringBuffer publString = new StringBuffer("JRNL        PUBL   ");
        StringBuffer refnString = new StringBuffer("JRNL        REFN                   ");
        StringBuffer pmidString = new StringBuffer("JRNL        PMID   ");
        StringBuffer doiString = new StringBuffer("JRNL        DOI    ");

        for (Author author : authorList) {
            authString.append(author).append(",");
        }
        jrnlString.append(authString.toString()).append(eol);
        titlString.append(title);
        jrnlString.append(titlString.toString()).append(eol);
        if (!editorList.isEmpty()) {
            for (Author editor : editorList) {
                editString.append(editor).append(",");
            }
            jrnlString.append(editString.toString()).append(eol);
        }
//        if (ref.equals("TO BE PUBLISHED")) {
            refString.append(ref);
        jrnlString.append(refString.toString()).append(eol);
//        }
        if (!publisher.equals("")) {
            publString.append(publisher);
            jrnlString.append(publString.toString()).append(eol);
        }
        if (!refn.equals("")) {
            refnString.append(refn);
            jrnlString.append(refnString.toString()).append(eol);
        }
        if (!pmid.equals("")) {
            pmidString.append(pmid);
            jrnlString.append(pmidString.toString()).append(eol);
        }
        if (!doi.equals("")) {
            doiString.append(doi);
            jrnlString.append(doiString.toString()).append(eol);
        }

        return jrnlString.toString();
    }
}
