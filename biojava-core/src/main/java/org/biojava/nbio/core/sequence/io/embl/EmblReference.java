package org.biojava.nbio.core.sequence.io.embl;

/**
 * This class contains the processed data of embl file that
 * contains the referenceNumber, referenceComment, referencePosition
 * referenceCrossReference, referenceGroup, referenceAuthor
 * referenceTitle, referenceLocation
 * @since 5.0.0
 * @author Noor Aldeen Al Mbaidin
 */
public class EmblReference {

    
    private String referenceNumber;
    private String referenceComment;
    private String referencePosition;
    private String referenceCrossReference;
    private String referenceGroup;
    private String referenceAuthor;
    private String referenceTitle;
    private String referenceLocation;

    /**
     * The RN (Reference Number) line gives a unique number to each reference
     * Citation within an entry. This number is used to designate the reference
     * in comments and in the feature table.
     * @return referenceNumber
     */
    public String getReferenceNumber() {
        return referenceNumber;
    }

    public void setReferenceNumber(String referenceNumber) {
        this.referenceNumber = referenceNumber;
    }

    /**
     * The RC (Reference Comment) linetype is an optional linetype which appears if
     * The reference has a comment.
     * @return String
     */
    public String getReferenceComment() {
        return referenceComment;
    }

    public void setReferenceComment(String referenceComment) {
        this.referenceComment = referenceComment;
    }

    /**
     * The RP (Reference Position) linetype is
     * an optional linetype which appears if
     * one or more contiguous base spans of
     * the presented sequence can be attributed
     * to the reference in question.
     * @return String
     */
    public String getReferencePosition() {
        return referencePosition;
    }

    public void setReferencePosition(String referencePosition) {
        this.referencePosition = referencePosition;
    }

    /**
     * The RX (reference cross-reference) linetype is
     * an optional linetype which appears if
     * one or more contiguous base spans of the
     * presented sequence can be attributed
     * to the reference in question.
     * @return String
     */
    public String getReferenceCrossReference() {
        return referenceCrossReference;
    }

    public void setReferenceCrossReference(String referenceCrossReference) {
        this.referenceCrossReference = referenceCrossReference;
    }

    /**
     * The RG (Reference Group) lines list the working groups/consortia that
     * produced the record.
     * @return String
     */
    public String getReferenceGroup() {
        return referenceGroup;
    }

    public void setReferenceGroup(String referenceGroup) {
        this.referenceGroup = referenceGroup;
    }

    /**
     * The RA (Reference Author) lines list the authors of the paper (or other
     * work) cited. All of the authors are included, and are listed in the order
     * given in the paper.
     * @return String
     */
    public String getReferenceAuthor() {
        return referenceAuthor;
    }

    public void setReferenceAuthor(String referenceAuthor) {
        this.referenceAuthor = referenceAuthor;
    }

    /**
     * The RT (Reference Title) lines give the title of the paper (or other work) as
     * exactly as is possible given the limitations of computer character sets.
     * @return String
     */
    public String getReferenceTitle() {
        return referenceTitle;
    }

    public void setReferenceTitle(String referenceTitle) {
        this.referenceTitle = referenceTitle;
    }

    /**
     * The RL (Reference Location) line contains the conventional citation
     * information for the reference.
     * @return String
     */
    public String getReferenceLocation() {
        return referenceLocation;
    }

    public void setReferenceLocation(String referenceLocation) {
        this.referenceLocation = referenceLocation;
    }
}
