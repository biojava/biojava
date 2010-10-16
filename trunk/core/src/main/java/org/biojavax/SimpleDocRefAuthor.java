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
 *      http://www.biojava.orDocRef
 */

package org.biojavax;
import org.biojava.utils.Unchangeable;


/**
 * Represents an author of a documentary reference.
 * @author Richard Holland
 * @author George Waldon
 * @see DocRef
 * @since 1.5
 */
public class SimpleDocRefAuthor extends Unchangeable implements DocRefAuthor  {
    
    private String name;
    private boolean editor;
    private boolean consortium;
    
    /**
     * Constructs a new author instance.
     * @param name the author name
     * @param consortium are they a consortium?
     * @param editor are they an editor?
     */
    public SimpleDocRefAuthor(String name, boolean consortium, boolean editor) {
        if (name==null) throw new IllegalArgumentException("Name cannot be null");
        this.name = name;
        this.editor = editor;
        if(editor==true) {
            int idx = name.lastIndexOf(" (ed.)");
            if(idx!=-1 && idx==(name.length()-6))
                this.name = name.substring(0,name.length()-6);
        }
        this.consortium = consortium;
        if(consortium==true) {
            int idx = name.lastIndexOf(" (consortium)");
            if(idx!=-1 && idx==(name.length()-13))
                this.name = name.substring(0,name.length()-13);
        }
    }
    
    /**
     * Constructs a new author instance from a string.
     * @param s the input string (author name with (ed.) and (consortium) suffixes).
     */
    public SimpleDocRefAuthor(String s) {
        if (s==null) throw new IllegalArgumentException("Name cannot be null");
        String[] parts = s.split("\\(");
        this.name = parts[0].trim();
        if (parts.length ==3) {
            parts[1] = parts[1].trim();
            parts[2] = parts[2].trim();
            parts[1] = parts[1].substring(0,parts[1].length()-1); // chomp bracket
            parts[2] = parts[2].substring(0,parts[2].length()-1); // chomp bracket
            if (parts[1].equals("ed.") || parts[2].equals("ed.")) this.editor=true;
            else this.editor = false;
            if (parts[1].equals("consortium") || parts[2].equals("consortium")) this.consortium = true;
            else this.consortium = false;
        } else if (parts.length ==2) {
            parts[1] = parts[1].trim();
            parts[1] = parts[1].substring(0,parts[1].length()-1); // chomp bracket
            if (parts[1].equals("ed.")) this.editor=true;
            else this.editor = false;
            if (parts[1].equals("consortium")) this.consortium = true;
            else this.consortium = false;
        } else {
            this.editor = false;
            this.consortium = false;
        }
        if(editor==true) {
            int idx = this.name.lastIndexOf(" (ed.)");
            if(idx!=-1 && idx==(this.name.length()-6))
                this.name = this.name.substring(0,this.name.length()-6);
        }
        if(consortium==true) {
            int idx = this.name.lastIndexOf(" (consortium)");
            if(idx!=-1 && idx==(this.name.length()-13))
                this.name = this.name.substring(0,this.name.length()-13);
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public String getName() {
        return this.name;
    }
    
    /**
     * {@inheritDoc}
     */
    public String getExtendedName() {
        StringBuffer result = new StringBuffer();
        result.append(this.name);
        if (this.consortium) result.append(" (consortium)");
        if (this.editor) result.append(" (ed.)");
        return result.toString();
    }
    
    /**
     * {@inheritDoc}
     */
    public boolean isEditor() {
        return this.editor;
    }
    
    /**
     * {@inheritDoc}
     */
    public boolean isConsortium() {
        return this.consortium;
    }
    
    /**
     * {@inheritDoc}
     * Document authors are compared first by name, then consortium status, then editor status.
     */
    public int compareTo(Object o) {
        if (o==this) return 0;
        DocRefAuthor them = (DocRefAuthor)o;
        if (!this.name.equals(them.getName())) return this.name.compareTo(them.getName());
        if (this.consortium!=them.isConsortium()) return this.consortium?-1:1;
        if (this.editor!=them.isEditor()) return this.editor?-1:1;
        return 0;
    }
    
    /**
     * {@inheritDoc}
     * Document references are equal if they have all fields the same.
     */
    public boolean equals(Object obj) {
        if(this == obj) return true;
        if (obj==null || !(obj instanceof DocRefAuthor)) return false;
        DocRefAuthor them = (DocRefAuthor)obj;
        return (this.name.equals(them.getName()) &&
                this.consortium==them.isConsortium() &&
                this.editor==them.isEditor());
    }
    
    /**
     * {@inheritDoc}
     */
    public int hashCode() {
        int code = 17;
        code = 37*code + this.name.hashCode();
        code = 37*code + (this.consortium?1:0);
        code = 37*code + (this.editor?1:0);
        return code;
    }
    
    /**
     * {@inheritDoc}
     * Form: "name (consortium) (ed.)" where sections in brackets are optional.
     */
    public String toString() {
        return this.getExtendedName();
    }
}