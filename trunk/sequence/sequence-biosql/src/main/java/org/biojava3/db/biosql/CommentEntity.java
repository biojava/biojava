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
package org.biojava3.db.biosql;

import java.io.Serializable;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.NamedQueries;
import javax.persistence.NamedQuery;
import javax.persistence.Table;

/**
 * Persitent entity bean for CommentEntity
 * @author Mark Schreiber
 */
@Entity
@Table(name = "COMMENT")
@NamedQueries({
    @NamedQuery(name = "CommentEntity.findByCommentId", query = "SELECT c FROM CommentEntity c WHERE c.commentId = :commentId"), 
    @NamedQuery(name = "CommentEntity.findByCommentText", query = "SELECT c FROM CommentEntity c WHERE c.commentText = :commentText"), 
    @NamedQuery(name = "CommentEntity.findByRank", query = "SELECT c FROM CommentEntity c WHERE c.rank = :rank")})
public class CommentEntity implements Serializable {
    private static final long serialVersionUID = 1L;
    @Id
    @GeneratedValue(strategy=GenerationType.IDENTITY)
    @Column(name = "COMMENT_ID", nullable = false)
    private Integer commentId;
    @Column(name = "COMMENT_TEXT", nullable = false)
    private String commentText;
    @Column(name = "RANK", nullable = false)
    private int rank;
    @JoinColumn(name = "BIOENTRY_ID", referencedColumnName = "BIOENTRY_ID")
    @ManyToOne
    private BioentryEntity bioentryId;

    public CommentEntity() {
    }

    public CommentEntity(Integer commentId) {
        this.commentId = commentId;
    }

    public CommentEntity(Integer commentId, String commentText, int rank) {
        this.commentId = commentId;
        this.commentText = commentText;
        this.rank = rank;
    }

    public Integer getCommentId() {
        return commentId;
    }

    public void setCommentId(Integer commentId) {
        this.commentId = commentId;
    }

    public String getCommentText() {
        return commentText;
    }

    public void setCommentText(String commentText) {
        this.commentText = commentText;
    }

    public int getRank() {
        return rank;
    }

    public void setRank(int rank) {
        this.rank = rank;
    }

    public BioentryEntity getBioentryId() {
        return bioentryId;
    }

    public void setBioentryId(BioentryEntity bioentryId) {
        this.bioentryId = bioentryId;
    }

    @Override
    public int hashCode() {
        int hash = 0;
        hash += (commentId != null ? commentId.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        // TODO: Warning - this method won't work in the case the id fields are not set
        if (!(object instanceof CommentEntity)) {
            return false;
        }
        CommentEntity other = (CommentEntity) object;
        if ((this.commentId == null && other.commentId != null) || (this.commentId != null && !this.commentId.equals(other.commentId))) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "org.biosql.entity.CommentEntity[commentId=" + commentId + "]";
    }

}
