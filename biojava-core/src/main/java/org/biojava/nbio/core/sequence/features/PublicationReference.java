package org.biojava.nbio.core.sequence.features;

import java.util.ArrayList;
import java.util.List;

public class PublicationReference {

    public enum ReferenceType {
       UNKNOWN, PUBMED, PATENT, DIRECT_SUBMISSION;
    }

    private String id, title, journal;
    private ReferenceType referenceType;
    private List<PublicationReferenceAuthor> authors = new ArrayList<>();

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getTitle() {
        return title;
    }

    public void setTitle(String title) {
        this.title = title;
    }

    public String getJournal() {
        return journal;
    }

    public void setJournal(String journal) {
        this.journal = journal;
    }

    public ReferenceType getReferenceType() {
        return referenceType;
    }

    public void setReferenceType(ReferenceType referenceType) {
        this.referenceType = referenceType;
    }

    public List<PublicationReferenceAuthor> getAuthors() {
        return authors;
    }

    public void setAuthors(List<PublicationReferenceAuthor> authors) {
        this.authors = authors;
    }

    @Override
    public String toString() {
        return "PublicationReference{" +
                "id='" + id + '\'' +
                ", title='" + title + '\'' +
                ", journal='" + journal + '\'' +
                ", referenceType=" + referenceType +
                ", authors=" + authors +
                '}';
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        PublicationReference that = (PublicationReference) o;

        if (id != null ? !id.equals(that.id) : that.id != null) return false;
        if (title != null ? !title.equals(that.title) : that.title != null) return false;
        if (journal != null ? !journal.equals(that.journal) : that.journal != null) return false;
        if (referenceType != that.referenceType) return false;
        return authors != null ? authors.equals(that.authors) : that.authors == null;

    }

    @Override
    public int hashCode() {
        int result = id != null ? id.hashCode() : 0;
        result = 31 * result + (title != null ? title.hashCode() : 0);
        result = 31 * result + (journal != null ? journal.hashCode() : 0);
        result = 31 * result + (referenceType != null ? referenceType.hashCode() : 0);
        result = 31 * result + (authors != null ? authors.hashCode() : 0);
        return result;
    }
}