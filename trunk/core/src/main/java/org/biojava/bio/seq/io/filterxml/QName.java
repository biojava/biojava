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

package org.biojava.bio.seq.io.filterxml;

class QName {
    String nsURI;
    String localName;
    
    public QName(String nsURI, String localName) {
        this.nsURI = nsURI;
        this.localName = localName;
    }
    
    public int hashCode() {
        return (nsURI == null ? 0 : nsURI.hashCode()) + localName.hashCode();
    }
    
    public boolean equals(Object o) {
        if (o instanceof QName) {
            QName qno = (QName) o;
            if (qno.nsURI == null) {
                if (nsURI != null) {
                    return false;
                }
            } else {
                if (!qno.nsURI.equals(nsURI)) {
                    return false;
                }
            }
            return qno.localName.equals(localName);
        }
        return false;
    }
}
